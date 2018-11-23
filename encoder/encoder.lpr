program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, extern, ap, conv, correlation, anysort;

const
  BandCount = 4;
  C1Freq = 32.703125;

type
  TEncoder = class;
  TFrame = class;
  TBand = class;
  TChunk = class;

  TBandGlobalData = record
    fcl, fch: Double;
    underSample: Integer;
    filteredData: TDoubleDynArray;
    dstData: TSmallIntDynArray;
  end;

  PBandGlobalData = ^TBandGlobalData;

  { TChunk }

  TChunkList = specialize TFPGObjectList<TChunk>;

  TChunk = class
  public
    frame: TFrame;
    reducedChunk: TChunk;

    index: Integer;
    underSample: Integer;
    dstBitShift: Integer;

    origSrcData: PDouble;
    srcData: TDoubleDynArray;
    dct: TDoubleDynArray;
    dstData: TByteDynArray;

    constructor Create(frm: TFrame; idx: Integer; underSmp: Integer; srcDta: PDouble);
    destructor Destroy; override;

    procedure ComputeDCT;
    procedure ComputeBitRange;
    procedure MakeSrcData(origData: PDouble);
    procedure MakeDstData;
  end;

  { TBand }

  TBand = class
  public
    frame: TFrame;

    index: Integer;
    ChunkCount: Integer;

    srcData: PDouble;
    dstData: TDoubleDynArray;

    finalChunks: TChunkList;

    globalData: PBandGlobalData;

    constructor Create(frm: TFrame; idx: Integer; startSample, endSample: Integer);
    destructor Destroy; override;

    procedure MakeChunks;
    procedure FindBestFinalChunks;
    procedure MakeDstData;
  end;

  { TFrame }

  TFrame = class
  public
    encoder: TEncoder;

    index: Integer;
    SampleCount: Integer;
    FrameSize: Integer;

    chunkRefs, reducedChunks: TChunkList;
    Centroids: TStringList;

    bands: array[0..BandCount - 1] of TBand;

    constructor Create(enc: TEncoder; idx, startSample, endSample: Integer);
    destructor Destroy; override;

    procedure MakeChunks;
    procedure KMeansReduce;
    procedure MakeDstData;
  end;

  TFrameList = specialize TFPGObjectList<TFrame>;

  { TEncoder }

  TEncoder = class
  public
    inputFN, outputFN: String;

    BitRate: Integer;
    Precision: Integer;
    BandTransFactor: Double;
    LowCut: Double;
    HighCut: Double;
    ChunkBitDepth: Integer; // 1 to 8 Bits
    ChunkSize: Integer;
    MaxChunksPerFrame: Integer;
    ReduceBassBand: Boolean;
    VariableFrameSizeRatio: Double;
    TrebleBoost: Boolean;

    SampleRate: Integer;
    SampleCount: Integer;
    BlockSampleCount: Integer;
    ProjectedByteSize, FrameCount: Integer;
    Verbose: Boolean;

    srcHeader: array[$00..$2b] of Byte;
    srcData: TSmallIntDynArray;
    dstData: TSmallIntDynArray;

    frames: TFrameList;

    bandData: array[0 .. BandCount - 1] of TBandGlobalData;

    class function make16BitSample(smp: Double): SmallInt;
    class function makeOutputSample(smp: Double; OutBitDepth, bitShift: Integer): Byte;
    class function makeFloatSample(smp: SmallInt): Double; overload;
    class function makeFloatSample(smp: Byte; OutBitDepth, bitShift: Integer): Double; overload;
    class function ComputeDCT(chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
    class function ComputeInvDCT(chunkSz: Integer; const dct: TDoubleDynArray): TDoubleDynArray;
    class function CompareDCT(firstCoeff, lastCoeff: Integer; const dctA, dctB: TDoubleDynArray): Double;
    class function CheckJoinPenalty(x, y, z, a, b, c: Double; TestRange: Boolean): Boolean; inline;

    constructor Create(InFN, OutFN: String);
    destructor Destroy; override;

    procedure Load;
    procedure SaveWAV;
    procedure SaveRSC;
    procedure SaveStream(AStream: TStream);
    procedure SaveBandWAV(index: Integer; fn: String);

    procedure MakeBandGlobalData;
    procedure MakeBandSrcData(AIndex: Integer);

    procedure PrepareFrames;
    procedure MakeFrames;
    procedure MakeDstData;

    function DoFilterCoeffs(fc, transFactor: Double; HighPass, Windowed: Boolean): TDoubleDynArray;
    function DoFilter(const samples, coeffs: TDoubleDynArray): TDoubleDynArray;
    function DoBPFilter(fcl, fch, transFactor: Double; const samples: TDoubleDynArray): TDoubleDynArray;

    function ComputeEAQUAL(chunkSz: Integer; UseDIX, Verbz: Boolean; const smpRef, smpTst: TSmallIntDynArray): Double;
    function ComputeCorrelation(const smpRef, smpTst: TSmallIntDynArray): Double;
  end;


function IsDebuggerPresent(): LongBool stdcall; external 'kernel32.dll';

function HasParam(p: String): Boolean;
var i: Integer;
begin
  Result := False;
  for i := 3 to ParamCount do
    if SameText(p, ParamStr(i)) then
      Exit(True);
end;

function ParamStart(p: String): Integer;
var i: Integer;
begin
  Result := -1;
  for i := 3 to ParamCount do
    if AnsiStartsStr(p, ParamStr(i)) then
      Exit(i);
end;

function ParamValue(p: String; def: Double): Double;
var
  idx: Integer;
begin
  idx := ParamStart(p);
  if idx < 0 then
    Exit(def);
  Result := StrToFloatDef(copy(ParamStr(idx), Length(p) + 1), def);
end;

{ TChunk }

constructor TChunk.Create(frm: TFrame; idx: Integer; underSmp: Integer; srcDta: PDouble);
begin
  index := idx;
  underSample := underSmp;
  frame := frm;

  SetLength(srcData, frame.encoder.chunkSize);

  if Assigned(srcDta) then
  begin
    origSrcData := @srcDta[idx * frame.encoder.chunkSize * underSample];
    MakeSrcData(origSrcData);
  end;

  reducedChunk := Self;
end;

destructor TChunk.Destroy;
begin
  inherited Destroy;
end;

procedure TChunk.ComputeDCT;
var
  i: Integer;
  data: TDoubleDynArray;
begin
  SetLength(data, Length(srcData));
  for i := 0 to High(data) do
    data[i] := srcData[i] * (1 shl dstBitShift);
  dct := TEncoder.ComputeDCT(Length(data), data);
end;

procedure TChunk.ComputeBitRange;
var
  i, hiSmp: Integer;
  bitRange: Integer;
begin
  hiSmp := 0;
  for i := 0 to frame.encoder.chunkSize - 1 do
    hiSmp := max(hiSmp, abs(TEncoder.make16BitSample(srcData[i])));
  bitRange := EnsureRange(1 + ceil(log2(hiSmp + 1.0)), 8, 16);
  dstBitShift := 16 - bitRange;
end;

procedure TChunk.MakeSrcData(origData: PDouble);
var
  j, k, pos, n: Integer;
  acc: Double;
begin
  for j := 0 to frame.encoder.chunkSize - 1 do
  begin
    pos := j * underSample;

    acc := 0.0;
    n := 0;
    for k := 0 to underSample - 1 do
    begin
      if pos + k >= frame.SampleCount then
        Break;
      acc += origData[pos + k];
      Inc(n);
    end;

    if n = 0 then
      srcData[j] := 0
    else
      srcData[j] := acc / n;
  end;
end;

procedure TChunk.MakeDstData;
var
  i: Integer;
begin
  SetLength(dstData, frame.encoder.chunkSize);
  for i := 0 to frame.encoder.chunkSize - 1 do
    dstData[i] := TEncoder.makeOutputSample(srcData[i], frame.encoder.ChunkBitDepth, dstBitShift);
end;

{ TBand }

constructor TBand.Create(frm: TFrame; idx: Integer; startSample, endSample: Integer);
begin
  frame := frm;
  index := idx;
  globalData := @frame.encoder.bandData[index];

  srcData := @globalData^.filteredData[startSample];

  ChunkCount := (endSample - startSample + 1 - 1) div (frame.encoder.ChunkSize * globalData^.underSample) + 1;

  finalChunks := TChunkList.Create;
end;

destructor TBand.Destroy;
begin
  finalChunks.Free;

  inherited Destroy;
end;

procedure TBand.MakeChunks;
var
  i: Integer;
  chunk: TChunk;
begin
  finalChunks.Clear;
  finalChunks.Capacity := ChunkCount;
  for i := 0 to ChunkCount - 1 do
  begin
    chunk := TChunk.Create(frame, i, globalData^.underSample, srcData);
    chunk.ComputeBitRange;
    chunk.ComputeDCT;
    chunk.MakeDstData;
    finalChunks.Add(chunk);
  end;
end;

procedure TBand.FindBestFinalChunks;
var
  i, j, prec: Integer;
  Clusters: TIntegerDynArray;
  Dataset: TReal2DArray;
begin
  prec := frame.encoder.Precision;
  if (prec = 0) or (not frame.encoder.ReduceBassBand and (index = 0)) then Exit;

  SetLength(Dataset, finalChunks.Count, frame.encoder.chunkSize);

  for i := 0 to finalChunks.Count - 1 do
    for j := 0 to frame.encoder.chunkSize - 1 do
      Dataset[i, j] := finalChunks[i].dct[j];

  Clusters := nil;
  DoExternalYakmo(Dataset, 0, 1, True, False, frame.Centroids, Clusters);

  for i := 0 to finalChunks.Count - 1 do
    finalChunks[i].reducedChunk := frame.reducedChunks[Clusters[i]];
end;

procedure TBand.MakeDstData;
var
  i, j, k, pos: Integer;
  chunk: TChunk;
  smp: Double;
begin
  //WriteLn('MakeDstData #', index);

  SetLength(dstData, frame.SampleCount);
  FillQWord(dstData[0], frame.SampleCount, 0);

  pos := 0;
  for i := 0 to finalChunks.Count - 1 do
  begin
    chunk := finalChunks[i];

    for j := 0 to frame.encoder.chunkSize - 1 do
    begin
      smp := TEncoder.makeFloatSample(chunk.reducedChunk.dstData[j], frame.encoder.ChunkBitDepth, chunk.dstBitShift);

      for k := 0 to globalData^.underSample - 1 do
      begin
        if InRange(pos, 0, High(dstData)) then
          dstData[pos] := smp;
        Inc(pos);
      end;
    end;
  end;
end;

constructor TFrame.Create(enc: TEncoder; idx, startSample, endSample: Integer);
var
  i: Integer;
begin
  encoder := enc;
  index := idx;
  SampleCount := endSample - startSample + 1;

  //Assert(ChunkCount < 256 * MDBlockBufferLen * trunc(8 / log2(encoder.ChunkBitDepth)), 'Frame too big! (VariableFrameSizeRatio too high and/or BitRate too low)');

  for i := 0 to BandCount - 1 do
    bands[i] := TBand.Create(Self, i, startSample, endSample);

  if encoder.Verbose then
  begin
    Write('Frame #', index);
    for i := 0 to BandCount - 1 do
      Write(#9, bands[i].ChunkCount);
    WriteLn;
  end;

  chunkRefs := TChunkList.Create(False);
  reducedChunks := TChunkList.Create;
  Centroids := TStringList.Create;
end;

destructor TFrame.Destroy;
var
  i: Integer;
begin
  Centroids.Free;
  reducedChunks.Free;
  chunkRefs.Free;

  for i := 0 to BandCount - 1 do
    bands[i].Free;

  inherited Destroy;
end;

procedure TFrame.MakeDstData;
var
  i: Integer;
begin
  for i := 0 to BandCount - 1 do
  begin
    bands[i].FindBestFinalChunks;
    bands[i].MakeDstData;
  end;
end;

procedure TFrame.MakeChunks;
var
  i, j, k: Integer;
begin
  chunkRefs.Clear;
  for i := 0 to BandCount - 1 do
  begin
    bands[i].MakeChunks;
    for j := Ord(not encoder.ReduceBassBand) to bands[i].finalChunks.Count - 1 do
      for k := 1 to round(power(bands[i].globalData^.underSample, sqrt(2.0))) do
        chunkRefs.Add(bands[i].finalChunks[j]);
  end;
end;

procedure TFrame.KMeansReduce;
var
  i, j, prec: Integer;
  chunk: TChunk;
  centroid: TDoubleDynArray;
  Clusters: TIntegerDynArray;
  Dataset: TReal2DArray;
begin
  prec := encoder.Precision;
  if (prec = 0) or (chunkRefs.Count <= encoder.MaxChunksPerFrame) then Exit;

  if encoder.Verbose then
    WriteLn('KMeansReduce Frame = ', index, ', N = ', chunkRefs.Count);

  SetLength(Dataset, chunkRefs.Count, encoder.chunkSize);

  for i := 0 to chunkRefs.Count - 1 do
    for j := 0 to encoder.chunkSize - 1 do
      Dataset[i, j] := chunkRefs[i].dct[j];

  Clusters := nil;
  DoExternalYakmo(Dataset, encoder.MaxChunksPerFrame, prec, False, False, Centroids, Clusters);

  reducedChunks.Clear;
  reducedChunks.Capacity := encoder.MaxChunksPerFrame;
  for i := 0 to encoder.MaxChunksPerFrame - 1 do
  begin
    chunk := TChunk.Create(Self, i, 1, nil);

    reducedChunks.Add(chunk);

    centroid := GetSVMLightLine(i, Centroids);

    centroid := TEncoder.ComputeInvDCT(encoder.ChunkSize, centroid);

    SetLength(chunk.dstData, encoder.chunkSize);
    for j := 0 to encoder.chunkSize - 1 do
      chunk.dstData[j] := TEncoder.makeOutputSample(centroid[j], encoder.ChunkBitDepth, 0);
  end;
end;

{ TEncoder }

procedure TEncoder.Load;
var
  fs: TFileStream;
  i: Integer;
begin
  WriteLn('Load ', inputFN);
  fs := TFileStream.Create(inputFN, fmOpenRead or fmShareDenyNone);
  try
    fs.ReadBuffer(srcHeader[0], SizeOf(srcHeader));
    SampleCount := (fs.Size - fs.Position) div SizeOf(SmallInt);
    SetLength(srcData, SampleCount);
    FillWord(srcData[0], SampleCount, 0);
    for i := 0 to SampleCount - 1 do
      srcData[i] := SmallInt(fs.ReadWord);
  finally
    fs.Free;
  end;

  SampleRate := PInteger(@srcHeader[$18])^;
end;

procedure TEncoder.SaveWAV;
var
  fs: TFileStream;
  wavFN: String;
begin
  wavFN := ChangeFileExt(outputFN, '.wav');

  WriteLn('Save ', wavFN);

  fs := TFileStream.Create(wavFN, fmCreate or fmShareDenyWrite);
  try
    fs.WriteBuffer(srcHeader[0], SizeOf(srcHeader));
    fs.WriteBuffer(dstData[0], SampleCount * 2);
  finally
    fs.Free;
  end;
end;

procedure TEncoder.SaveRSC;
var
  fs: TFileStream;
  cur: TMemoryStream;
  fn: String;
begin
  fs := nil;
  fn := ChangeFileExt(outputFN, '.rsc');
  cur := TMemoryStream.Create;
  fs := TFileStream.Create(fn, fmCreate);
  try
    WriteLn('Save ', fn);

    SaveStream(cur);
    cur.Position := 0;

    fs.CopyFrom(cur, cur.Size);

    writeln('FinalByteSize = ', cur.Size);
    writeln('FinalBitRate = ', round(cur.size * (8 / 1024) / (SampleCount / SampleRate)));
  finally
    fs.Free;
    cur.Free;
  end;
end;

procedure TEncoder.SaveStream(AStream: TStream);
//var
//  i, j, k, l: Integer;
//  a, b : Integer;
//  cl: TChunkList;
begin
  //for i := 0 to FrameCount - 1 do
  //begin
  //  for j := 0 to BandCount - 1 do
  //  begin
  //    Assert(frames[i].bands[j].chunks.Count mod (BlockSampleCount div ChunkSize) = 0);
  //    AStream.WriteByte(frames[i].bands[j].chunks.Count div (BlockSampleCount div ChunkSize));
  //
  //    cl := frames[i].bands[j].reducedChunks;
  //    if cl.Count = 0 then
  //      cl := frames[i].bands[j].chunks;
  //
  //    for l := 0 to ChunkSize - 1 do
  //      for k := 0 to cl.Count - 1 do
  //      begin
  //        if l < Length(cl[k].dstData) then
  //          AStream.WriteByte(cl[k].dstData[l])
  //        else
  //          AStream.WriteByte($80);
  //      end;
  //  end;
  //
  //  for j := 0 to BandCount - 1 do
  //  begin
  //    cl := frames[i].bands[j].chunks;
  //
  //    for k := 0 to cl.Count - 1 do
  //    begin
  //      if (k and 7 = 0) and (ChunkBitDepth > 8) then
  //      begin
  //        b := 0;
  //        for l := k to k + 7 do
  //        begin
  //          b := b shl 1;
  //          if l < cl.Count then
  //            b := b or (cl[l].dstBitShift shr 3);
  //        end;
  //        AStream.WriteByte(b);
  //      end;
  //
  //      a := cl[k].reducedChunk.index;
  //      AStream.WriteByte(a and $ff);
  //    end;
  //  end;
  //end;
  //
  //AStream.WriteByte(0); // Termination
end;

procedure TEncoder.SaveBandWAV(index: Integer; fn: String);
var
  fs: TFileStream;
begin
  //WriteLn('SaveBandWAV #', index, ' ', fn);

  fs := TFileStream.Create(fn, fmCreate or fmShareDenyWrite);
  try
    fs.WriteBuffer(srcHeader[0], SizeOf(srcHeader));
    fs.WriteBuffer(bandData[index].dstData[0], SampleCount * 2);
  finally
    fs.Free;
  end;
end;

procedure TEncoder.MakeBandGlobalData;
var
  i: Integer;
  ratio, hc: Double;
  bnd: TBandGlobalData;
begin
  for i := 0 to BandCount - 1 do
  begin
    bnd.dstData := nil;
    FillChar(bnd, SizeOf(bnd), 0);

    // determing low and high bandpass frequencies

    hc := min(HighCut, SampleRate / 2);
    ratio := (log2(hc) - log2(max(C1Freq, LowCut))) / BandCount;

    if i = 0 then
      bnd.fcl := LowCut / SampleRate
    else
      bnd.fcl := 0.5 * power(2.0, -floor((BandCount - i) * ratio));

    if i = BandCount - 1 then
      bnd.fch := hc / SampleRate
    else
      bnd.fch := 0.5 * power(2.0, -floor((BandCount - 1 - i) * ratio));

    // undersample if the band high freq is a lot lower than nyquist

    bnd.underSample := Max(1, round(0.25 / bnd.fch));

    bandData[i] := bnd;
  end;
end;

procedure TEncoder.MakeBandSrcData(AIndex: Integer);
var
  j: Integer;
  bnd: TBandGlobalData;
begin
  bnd := bandData[AIndex];

  SetLength(bnd.filteredData, SampleCount);
  for j := 0 to SampleCount - 1 do
    bnd.filteredData[j] := makeFloatSample(srcData[j]);

  // band pass the samples
  bnd.filteredData := DoBPFilter(bnd.fcl, bnd.fch, BandTransFactor, bnd.filteredData);
  SetLength(bnd.filteredData, SampleCount);

  bandData[AIndex] := bnd;
end;

procedure TEncoder.PrepareFrames;

  procedure DoBand(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  begin
    MakeBandSrcData(AIndex);
  end;

var
  i, k, fixedCost, frameCost, nextStart, psc: Integer;
  frm: TFrame;
  avgPower, totalPower, perFramePower, curPower: Double;
begin
  MakeBandGlobalData;

  // pass 1

  BlockSampleCount := 0;
  for i := 0 to BandCount - 1 do
    if bandData[i].underSample * ChunkSize > BlockSampleCount then
      BlockSampleCount := bandData[i].underSample * ChunkSize;

  // ensure srcData ends on a full block
  psc := SampleCount;
  SampleCount := ((SampleCount - 1) div BlockSampleCount + 1) * BlockSampleCount;
  SetLength(srcData, SampleCount);
  for i :=psc to SampleCount - 1 do
    srcData[i] := 0;

  ProjectedByteSize := ceil((SampleCount / SampleRate) * BitRate * (1024 / 8));

  frameCost := MaxChunksPerFrame * ChunkSize;
  fixedCost := 0;
  for i := 0 to BandCount - 1 do
    fixedCost += (SampleCount * (round(log2(MaxChunksPerFrame)) + round(log2(ChunkBitDepth)))) div (8 * ChunkSize * bandData[i].underSample) + SizeOf(Byte);

  FrameCount := (ProjectedByteSize - fixedCost - 1) div frameCost + 1;

  ProjectedByteSize := fixedCost + FrameCount * frameCost;

  writeln('SampleRate = ', SampleRate);
  writeln('FrameCount = ', FrameCount);

  Assert(FrameCount > 0, 'Negative FrameCount! (BitRate too low)');

  if Verbose then
  begin
    WriteLn('SampleCount = ', SampleCount);
    writeln('FrameSize = ', ProjectedByteSize div FrameCount);
    writeln('ProjectedByteSize = ', ProjectedByteSize);
    writeln('ChunkSize = ', ChunkSize);
  end;

  ProcThreadPool.DoParallelLocalProc(@DoBand, 0, BandCount - 1, nil);

  // pass 2

  avgPower := 0.0;
  for i := 0 to SampleCount - 1 do
    avgPower += abs(makeFloatSample(srcData[i]));
  avgPower /= SampleCount;

  totalPower := 0.0;
  for i := 0 to SampleCount - 1 do
    totalPower += avgPower - avgPower * VariableFrameSizeRatio + abs(VariableFrameSizeRatio * makeFloatSample(srcData[i]));

  perFramePower := totalPower / FrameCount;

  if Verbose then
  begin
    writeln('TotalPower = ', FormatFloat('0.00', totalPower));
    writeln('PerFramePower = ', FormatFloat('0.00', perFramePower));
  end;

  k := 0;
  nextStart := 0;
  curPower := 0.0;
  for i := 0 to SampleCount - 1 do
  begin
    curPower += avgPower - avgPower * VariableFrameSizeRatio + abs(VariableFrameSizeRatio * makeFloatSample(srcData[i]));

    if (i mod BlockSampleCount = 0) and (curPower >= perFramePower) then
    begin
      frm := TFrame.Create(Self, k, nextStart, i - 1);
      frames.Add(frm);

      curPower := 0.0;
      nextStart := i;
      Inc(k);
    end;
  end;

  frm := TFrame.Create(Self, k, nextStart, SampleCount - 1);
  frames.Add(frm);

  FrameCount := frames.Count;

  for i := 0 to BandCount - 1 do
     WriteLn('Band #', i, ' (', round(bandData[i].fcl * SampleRate), ' Hz .. ', round(bandData[i].fch * SampleRate), ' Hz); ', bandData[i].underSample);
end;

procedure TEncoder.MakeFrames;

  procedure DoFrame(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    frm: TFrame;
  begin
    frm := frames[AIndex];

    frm.MakeChunks;
    frm.KMeansReduce;
    //frm.SortAndReindexReducedChunks; //TODO: maybe useful for better LZ compression
    frm.MakeDstData;
    Write('.');
  end;
begin
  ProcThreadPool.DoParallelLocalProc(@DoFrame, 0, FrameCount - 1, nil);
  WriteLn;
end;

function TEncoder.DoFilter(const samples, coeffs: TDoubleDynArray): TDoubleDynArray;
var
  i: Integer;
begin
  Result := nil;
  ConvR1D(samples, Length(samples), coeffs, Length(coeffs), Result);

  for i := 0 to High(samples) do
    Result[i] := Result[i + High(coeffs) div 2];

  SetLength(Result, Length(samples));
end;

function TEncoder.DoBPFilter(fcl, fch, transFactor: Double; const samples: TDoubleDynArray
  ): TDoubleDynArray;
var
  coeffs: TDoubleDynArray;
begin
  Result := samples;

  if fcl > 0.0 then
  begin
    coeffs := DoFilterCoeffs(fcl, transFactor, True, True);
    Result := DoFilter(Result, coeffs);
  end;

  if fch < 0.5 then
  begin
    coeffs := DoFilterCoeffs(fch, transFactor, False, True);
    Result := DoFilter(Result, coeffs);
  end;
end;

constructor TEncoder.Create(InFN, OutFN: String);
begin
  inputFN := InFN;
  outputFN := OutFN;

  BitRate := 128;
  Precision := 2;
  LowCut := 0.0;
  HighCut := 24000.0;
  ChunkBitDepth := 8;
  ChunkSize := 4;
  ReduceBassBand := False;
  TrebleBoost := False;
  VariableFrameSizeRatio := 0.0;

  MaxChunksPerFrame := 4096;
  BandTransFactor := 0.25;

  frames := TFrameList.Create;
end;

destructor TEncoder.Destroy;
begin
  frames.Free;

  inherited Destroy;
end;

procedure TEncoder.MakeDstData;
var
  i, j, k, pos: Integer;
  smp: Double;
  bnd: TBandGlobalData;
  resamp: array[0 .. BandCount - 1] of TDoubleDynArray;
begin
  WriteLn('MakeDstData');

  SetLength(dstData, SampleCount);
  FillWord(dstData[0], Length(dstData), 0);

  for i := 0 to BandCount - 1 do
  begin
    bnd := bandData[i];

    SetLength(bnd.dstData, Length(dstData));
    FillWord(bnd.dstData[0], Length(dstData), 0);

    bandData[i] := bnd;
  end;

  pos := 0;
  for k := 0 to frames.Count - 1 do
  begin
    for j := 0 to BandCount - 1 do
    begin
      bnd := bandData[j];
      resamp[j] := DoBPFilter(bnd.fcl, bnd.fch, BandTransFactor, frames[k].bands[j].dstData);
    end;

    for i := 0 to frames[k].SampleCount - 1 do
    begin
      for j := 0 to BandCount - 1 do
      begin
        smp := resamp[j][i];

        if InRange(pos, 0, High(dstData)) then
        begin
          bandData[j].dstData[pos] := make16BitSample(smp);
          dstData[pos] := EnsureRange(dstData[pos] + bandData[j].dstData[pos], Low(SmallInt), High(SmallInt));
        end;
      end;

      Inc(pos);
    end;
  end;
end;

function TEncoder.DoFilterCoeffs(fc, transFactor: Double; HighPass, Windowed: Boolean): TDoubleDynArray;
var
  sinc, win, sum: Double;
  i, N: Integer;
begin
  N := ceil(4.6 / (transFactor * fc));
  if (N mod 2) = 0 then N += 1;

  //writeln('DoFilterCoeffs ', ifthen(HighPass, 'HP', 'LP'), ' ', FloatToStr(SampleRate * fc), ' ', N);

  SetLength(Result, N);
  sum := 0;
  for i := 0 to N - 1 do
  begin
    sinc := 2.0 * fc * (i - (N - 1) / 2.0) * pi;
    if sinc = 0 then
      sinc := 1.0
    else
      sinc := sin(sinc) / sinc;

    win := 1.0;
    if Windowed then
    begin
{$if true}
      // blackman window
      win := 7938/18608 - 9240/18608 * cos(2 * pi * i / (N - 1)) + 1430/18608 * cos(4 * pi * i / (N - 1));
{$else}
      // sinc window
      win := (2 * i / (N - 1) - 1) * pi;
      if win = 0 then
        win := 1.0
      else
        win := sin(win) / win;
{$endif}
    end;

    Result[i] := sinc * win;
    sum += Result[i];
  end;

  if HighPass then
  begin
    for i := 0 to N - 1 do
      Result[i] := -Result[i] / sum;

    Result[(N - 1) div 2] += 1.0;
  end
  else
  begin
    for i := 0 to N - 1 do
      Result[i] := Result[i] / sum;
  end;
end;

class function TEncoder.make16BitSample(smp: Double): SmallInt;
begin
  Result := EnsureRange(round(smp * -Low(SmallInt)), Low(SmallInt), High(SmallInt));
end;

class function TEncoder.makeOutputSample(smp: Double; OutBitDepth, bitShift: Integer): Byte;
var
  smpI: Integer;
begin
  smpI := round(smp * (1 shl (OutBitDepth - 1 + bitShift)));
  Result := EnsureRange(smpI + (1 shl (OutBitDepth - 1)), 0, (1 shl OutBitDepth) - 1);
end;

class function TEncoder.makeFloatSample(smp: SmallInt): Double;
begin
  Result := smp / -Low(SmallInt);
end;

class function TEncoder.makeFloatSample(smp: Byte; OutBitDepth, bitShift: Integer): Double;
var
  smpI: Integer;
begin
  smpI := Integer(smp) - (1 shl (OutBitDepth - 1));
  Result := smpI / (1 shl (OutBitDepth - 1 + bitShift));
end;

class function TEncoder.ComputeDCT(chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
var
  k, n: Integer;
  sum, s: Double;
begin
  SetLength(Result, length(samples));
  for k := 0 to chunkSz - 1 do
  begin
    sum := 0;
    s := ifthen(k = 0, sqrt(0.5), 1.0);
    for n := 0 to chunkSz - 1 do
      sum += s * samples[n] * cos(pi * (n + 0.5) * k / chunkSz);
    Result[k] := sum * sqrt (2.0 / chunkSz);
  end;
end;

class function TEncoder.ComputeInvDCT(chunkSz: Integer; const dct: TDoubleDynArray): TDoubleDynArray;
var
  k, n: Integer;
  sum, s: Double;
begin
  SetLength(Result, length(dct));
  for n := 0 to chunkSz - 1 do
  begin
    sum := 0;
    for k := 0 to chunkSz - 1 do
    begin
      s := ifthen(k = 0, sqrt(0.5), 1.0);
      sum += s * dct[k] * cos (pi * (n + 0.5) * k / chunkSz);
    end;
    Result[n] := sum * sqrt(2.0 / chunkSz);
  end;
end;

class function TEncoder.CompareDCT(firstCoeff, lastCoeff: Integer; const dctA, dctB: TDoubleDynArray): Double;
var
  i: Integer;
begin
  Result := 0.0;

  for i := firstCoeff to lastCoeff do
    Result += sqr(dctA[i] - dctB[i]);

  Result := sqrt(Result);
end;

class function TEncoder.CheckJoinPenalty(x, y, z, a, b, c: Double; TestRange: Boolean): Boolean;
var
  dStart, dEnd: Double;
begin
  dStart := -1.5 * x + 2.0 * y - 0.5 * z;
  dEnd := -1.5 * a + 2.0 * b - 0.5 * c;

  Result := Sign(dStart) * Sign(dEnd) <> -1;
  if TestRange and Result then
    Result := InRange(y, a, c) or InRange(y, c, a);
end;

function TEncoder.ComputeEAQUAL(chunkSz: Integer; UseDIX, Verbz: Boolean; const smpRef, smpTst: TSmallIntDynArray): Double;
var
  i, DupCnt: Integer;
  FNRef, FNTst: String;
  ms: TMemoryStream;
begin
  DupCnt := Max(1, round(15 * SampleRate / chunkSz));

  FNRef := GetTempFileName('', 'ref-'+IntToStr(GetCurrentThreadId)+'.wav');
  FNTst := GetTempFileName('', 'tst-'+IntToStr(GetCurrentThreadId)+'.wav');

  ms := TMemoryStream.Create;
  try
    ms.Write(srcHeader[0], $28);
    ms.WriteDWord(chunkSz * SizeOf(SmallInt) * DupCnt);
    for i := 0 to DupCnt - 1  do
      ms.Write(smpRef[0], chunkSz * SizeOf(SmallInt));

    if SampleRate < 44100 then
    begin
      ms.Position := $18;
      ms.WriteDWord(44100);
    end;

    ms.SaveToFile(FNRef);
    ms.Clear;

    ms.Write(srcHeader[0], $28);
    ms.WriteDWord(chunkSz * SizeOf(SmallInt) * DupCnt);
    for i := 0 to DupCnt - 1  do
      ms.Write(smpTst[0], chunkSz * SizeOf(SmallInt));

    if SampleRate < 44100 then
    begin
      ms.Position := $18;
      ms.WriteDWord(44100);
    end;

    ms.SaveToFile(FNTst);
  finally
    ms.Free;
  end;

  Result := DoExternalEAQUAL(FNRef, FNTst, Verbz, UseDIX, -1);

  DeleteFile(FNRef);
  DeleteFile(FNTst);
end;

function TEncoder.ComputeCorrelation(const smpRef, smpTst: TSmallIntDynArray
  ): Double;
var
  i, len: Integer;
  rr, rt: TReal1DArray;
begin
  len := length(smpRef);
  Assert(len = length(smpTst), 'ComputeCorrelation length mismatch!');
  SetLength(rr, len);
  SetLength(rt, len);
  for i := 0 to len - 1 do
  begin
    rr[i] := makeFloatSample(smpRef[i]);
    rt[i] := makeFloatSample(smpTst[i]);
  end;

  Result := SpearmanRankCorrelation(rr, rt, len);
end;

procedure test_makeSample;
var
  i: Integer;
  b, bb, bbb, o: byte;
  f: Double;
begin
  for i := 0 to 255 do
  begin
    b := i mod 256;//random(256);
    bb := 8;//RandomRange(1, 9);
    bbb := 0;//RandomRange(0, 8);
    b := b and ((1 shl bb) - 1);

    f := TEncoder.makeFloatSample(b, bb, bbb);
    o := TEncoder.makeOutputSample(f, bb, bbb);
    writeln(b,#9,bb,#9,bbb,#9,o,#9,FloatToStr(f));
    assert(b = o);
  end;
end;

var
  enc: TEncoder;
  i: Integer;
  dix, cor: double;
  s: String;
begin
  try
    FormatSettings.DecimalSeparator := '.';

{$ifdef DEBUG}
    ProcThreadPool.MaxThreadCount := 1;
{$else}
    SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS);
{$endif}

    //test_makeSample;


    if ParamCount < 2 then
    begin
      WriteLn('Usage: ', ExtractFileName(ParamStr(0)) + ' <source file> <dest file> [options]');
      Writeln('Main options:');
      WriteLn(#9'-br'#9'encoder bit rate in kilobits/second; example: "-br128"');
      WriteLn(#9'-lc'#9'bass cutoff frequency');
      WriteLn(#9'-hc'#9'treble cutoff frequency');
      WriteLn(#9'-vfr'#9'RMS power based variable frame size ratio (0.0-1.0)');
      WriteLn(#9'-v'#9'verbose mode');
      Writeln('Development options:');
      WriteLn(#9'-cs'#9'chunk size');
      WriteLn(#9'-rbb'#9'enable lossy compression on bass band');
      WriteLn(#9'-cbd'#9'chunk bit depth (1-8)');
      WriteLn(#9'-pr'#9'K-means precision; 0: "lossless" mode');

      WriteLn;
      Writeln('(source file must be 16bit mono WAV)');
      WriteLn;
      Exit;
    end;

    enc := TEncoder.Create(ParamStr(1), ParamStr(2));
    try
      enc.BitRate := round(ParamValue('-br', enc.BitRate));
      enc.Precision := round(ParamValue('-pr', enc.Precision));
      enc.LowCut := ParamValue('-lc', enc.LowCut);
      enc.HighCut := ParamValue('-hc', enc.HighCut);
      enc.VariableFrameSizeRatio :=  EnsureRange(ParamValue('-vfr', enc.VariableFrameSizeRatio), 0.0, 1.0);
      enc.ChunkBitDepth := EnsureRange(round(ParamValue('-cbd', enc.ChunkBitDepth)), 1, 8);
      enc.ChunkSize := round(ParamValue('-cs', enc.ChunkSize));
      enc.Verbose := HasParam('-v');
      enc.ReduceBassBand := HasParam('-rbb');

      WriteLn('BitRate = ', FloatToStr(enc.BitRate));
      WriteLn('LowCut = ', FloatToStr(enc.LowCut));
      WriteLn('HighCut = ', FloatToStr(enc.HighCut));
      WriteLn('VariableFrameSizeRatio = ', FloatToStr(enc.VariableFrameSizeRatio));
      if enc.Verbose then
      begin
        WriteLn('ChunkSize = ', enc.ChunkSize);
        WriteLn('ReduceBassBand = ', BoolToStr(enc.ReduceBassBand, True));
        WriteLn('ChunkBitDepth = ', enc.ChunkBitDepth);
        WriteLn('Precision = ', enc.Precision);
      end;
      WriteLn;

      enc.Load;

      enc.PrepareFrames;
      enc.MakeFrames;
      enc.MakeDstData;

      enc.SaveWAV;
      if BandCount > 1 then
        for i := 0 to BandCount - 1 do
          enc.SaveBandWAV(i, ChangeFileExt(enc.outputFN, '-' + IntToStr(i) + '.wav'));
      if enc.Precision > 0 then
        enc.SaveRSC;

      dix := enc.ComputeEAQUAL(enc.SampleCount, False, True, enc.srcData, enc.dstData);
      cor := enc.ComputeCorrelation(enc.srcData, enc.dstData);

      WriteLn('EAQUAL = ', FloatToStr(dix));
      WriteLn('Correlation = ', FormatFloat('0.000000', cor));

      s := FloatToStr(dix) + ' ' + FormatFloat('0.000000', cor) + ' ';
      for i := 0 to ParamCount do s := s + ParamStr(i) + ' ';
      ShellExecute(0, 'open', 'cmd.exe', PChar('/c echo ' + s + ' >> ..\log.txt'), '', 0);

    finally
      enc.Free;
    end;

    WriteLn('Done.');
    if IsDebuggerPresent then
      ReadLn;

  except
    on e: Exception do
    begin
      WriteLn('Exception: ', e.Message, ' (', e.ClassName, ')');
      ReadLn;
    end;
  end;
end.

