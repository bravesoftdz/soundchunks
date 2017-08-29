program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, extern, ap, conv, anysort;

const
  BandCount = 1;
  MDSampleRate = 26390;

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

  { TCICFilter }

  TCICFilter = class
  private
    R, pos: Integer;
    integrator: TDoubleDynArray;
    comb: TDoubleDynArray2;
  public
    Ratio: Integer;
    Stages: Integer;
    HighPass: Boolean;

    constructor Create(AFc: Double; AStages: Integer; AHighPass: Boolean);

    function ProcessSample(Smp: Double): Double;
  end;

  { TChunk }

  TChunkList = specialize TFPGObjectList<TChunk>;

  TChunk = class
  public
    band: TBand;
    reducedChunk, sameChunk: TChunk;

    index: Integer;
    bitRange: Integer;
    dstBitShift: Integer;

    mixList: TChunkList;

    origSrcData: PDouble;
    srcData: TDoubleDynArray;
    dct: TDoubleDynArray;
    dstData: TByteDynArray;

    constructor Create(bnd: TBand; idx: Integer);
    destructor Destroy; override;

    procedure ComputeDCT;
    procedure InitMix;
    procedure AddToMix(chk: TChunk);
    procedure FinalizeMix;
    procedure ComputeBitRange;
    procedure MakeSrcData(origData: PDouble);
    procedure MakeDstData;
  end;

  { TBand }

  TBand = class
  public
    frame: TFrame;

    index: Integer;
    weight: Double;
    AlternateReduceMethod: Boolean;

    srcData: PDouble;
    dstData: TDoubleDynArray;

    chunks: TChunkList;
    reducedChunks: TChunkList;
    globalData: PBandGlobalData;

    constructor Create(frm: TFrame; idx: Integer; startSample, endSample: Integer);
    destructor Destroy; override;

    procedure MakeChunks;

    procedure KMeansReduce;
    procedure SortAndReindexReducedChunks;
    procedure MakeDstData;
  end;

  { TFrame }

  TFrame = class
  public
    encoder: TEncoder;

    index: Integer;
    ChunkCount: Integer;
    SampleCount: Integer;
    FrameSize: Integer;

    bands: array[0..BandCount - 1] of TBand;

    constructor Create(enc: TEncoder; idx, startSample, endSample: Integer);
    destructor Destroy; override;

    procedure MakeBands;
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
    StoredBitDepth: Integer; // max 8Bits
    OutputBitDepth: Integer; // 8 to 16Bits
    ChunkSize: Integer;
    ChunksPerBand: Integer;
    AlternateReduce: Boolean;
    BandBoundsOffset: Integer;
    BandDenoiseStages: Integer;
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
    class function CompareDCT(firstCoeff, lastCoeff: Integer; compress: Boolean; const dctA, dctB: TDoubleDynArray): Double;
    class function CompressDCT(coeff: Double): Double;
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
  end;


function IsDebuggerPresent () : LongBool stdcall; external 'kernel32.dll';

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

constructor TFrame.Create(enc: TEncoder; idx, startSample, endSample: Integer);
var
  i: Integer;
begin
  encoder := enc;
  index := idx;
  ChunkCount := (endSample - startSample + 1 - 1) div encoder.ChunkSize + 1;
  SampleCount := endSample - startSample + 1;

  Assert(ChunkCount < 65536 * trunc(8 / log2(encoder.OutputBitDepth - 7)), 'Frame too big! (VariableFrameSizeRatio too high and/or BitRate too low)');

  if encoder.Verbose then
  begin
    WriteLn('Frame #', index, #9, ChunkCount);
  end;

  for i := 0 to BandCount - 1 do
    bands[i] := TBand.Create(Self, i, startSample, endSample);
end;

destructor TFrame.Destroy;
var
  i: Integer;
begin
  for i := 0 to BandCount - 1 do
    bands[i].Free;

  inherited Destroy;
end;

procedure TFrame.MakeBands;
var
  i: Integer;
begin
  for i := 0 to BandCount - 1 do
    bands[i].MakeChunks;

  for i := 0 to BandCount - 1 do
  begin
    bands[i].KMeansReduce;
    bands[i].SortAndReindexReducedChunks;
    bands[i].MakeDstData;
  end;

  Write('.');
end;

constructor TCICFilter.Create(AFc: Double; AStages: Integer; AHighPass: Boolean);
begin
  HighPass := AHighPass;
  Stages := AStages;
  Ratio := round(0.5 / AFc);
  R := max(1, Ratio);
  SetLength(comb, Stages, R);
  SetLength(integrator, Stages);
  pos := 0;
end;

function TCICFilter.ProcessSample(Smp: Double): Double;
var
  i: Integer;
  lsmp: Double;
begin
  Result := 0;

  lsmp := comb[0, pos];

  for i := 0 to Stages - 1 do
  begin
    Result := Smp - comb[i, pos];
    comb[i, pos] := Smp;
    Smp := Result;
  end;

  pos := (pos + 1) mod R;

  for i := 0 to Stages - 1 do
  begin
    integrator[i] := integrator[i] + Result;
    Result := integrator[i];
  end;

  Result /= intpower(R, Stages);

  if HighPass then
    Result := lsmp - Result;
end;

constructor TBand.Create(frm: TFrame; idx: Integer; startSample, endSample: Integer);
begin
  frame := frm;
  index := idx;
  globalData := @frame.encoder.bandData[index];
  AlternateReduceMethod := frame.encoder.AlternateReduce;

  srcData := @globalData^.filteredData[startSample];

  chunks := TChunkList.Create;
  reducedChunks := TChunkList.Create;
end;

destructor TBand.Destroy;
begin
  chunks.Free;
  reducedChunks.Free;

  inherited Destroy;
end;

procedure TBand.MakeChunks;
var
  i: Integer;
  chunk: TChunk;
begin
  chunks.Clear;
  chunks.Capacity := (frame.chunkCount - 1) div globalData^.underSample + 1;
  for i := 0 to chunks.Capacity - 1 do
  begin
    chunk := TChunk.Create(Self, i);
    chunk.ComputeBitRange;
    chunk.MakeDstData;
    chunk.ComputeDCT;
    chunks.Add(chunk);
  end;
end;

procedure TBand.KMeansReduce;
var
  XYC: TIntegerDynArray;

  procedure DoXYC(AIndex: Integer);
  var
    i: Integer;
    reducedChunk: TChunk;
  begin
    reducedChunk := reducedChunks[AIndex];
    reducedChunk.InitMix;

    for i := 0 to chunks.Count - 1 do
      if XYC[i] = AIndex then
      begin
        reducedChunk.AddToMix(chunks[i]);
        chunks[i].reducedChunk := reducedChunk;
      end;

    reducedChunk.FinalizeMix;
  end;

var
  i, j, prec: Integer;
  chunk: TChunk;
  XY: TReal2DArray;
begin
  prec := frame.encoder.Precision;
  if (prec = 0) or (chunks.Count <= frame.encoder.ChunksPerBand) then Exit;

  //WriteLn('KMeansReduce #', index, ' ', globalData^.desiredChunkCount);

  SetLength(XY, chunks.Count, frame.encoder.chunkSize);

  for i := 0 to chunks.Count - 1 do
    for j := 0 to frame.encoder.chunkSize - 1 do
      XY[i, j] := chunks[i].dct[j];

  if AlternateReduceMethod then
    DoExternalSKLearn(XY, frame.encoder.ChunksPerBand, prec, False, XYC)
  else
    DoExternalYakmo(XY, frame.encoder.ChunksPerBand, 1, False, XYC);

  reducedChunks.Clear;
  reducedChunks.Capacity := frame.encoder.ChunksPerBand;
  for i := 0 to frame.encoder.ChunksPerBand - 1 do
  begin
    chunk := TChunk.Create(Self, i);

    reducedChunks.Add(chunk);
    DoXYC(i);
  end;
end;

function CompareReducedChunks(const Item1, Item2: TChunk): Integer;
begin
  Result := CompareValue(Item2.mixList.Count, Item1.mixList.Count);
end;

procedure TBand.SortAndReindexReducedChunks;
var i: Integer;
begin
  reducedChunks.Sort(@CompareReducedChunks);
  for i := 0 to reducedChunks.Count - 1 do
    reducedChunks[i].index := i;
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
  for i := 0 to chunks.Count - 1 do
  begin
    chunk := chunks[i];

    for j := 0 to frame.encoder.chunkSize - 1 do
    begin
      smp := TEncoder.makeFloatSample(chunk.reducedChunk.dstData[j], frame.encoder.StoredBitDepth, chunk.dstBitShift);

      for k := 0 to globalData^.underSample - 1 do
      begin
        if InRange(pos, 0, High(dstData)) then
          dstData[pos] := smp;
        Inc(pos);
      end;
    end;
  end;
end;

{ TChunk }

constructor TChunk.Create(bnd: TBand; idx: Integer);
begin
  index := idx;
  band := bnd;

  SetLength(srcData, band.frame.encoder.chunkSize);

  origSrcData := @band.srcData[idx * band.frame.encoder.chunkSize * band.globalData^.underSample];

  MakeSrcData(origSrcData);

  reducedChunk := Self;

  mixList := TChunkList.Create(False);
end;

destructor TChunk.Destroy;
begin
  mixList.Free;

  inherited Destroy;
end;

procedure TChunk.ComputeDCT;
var
  i: Integer;
  data: TDoubleDynArray;
begin
  SetLength(data, Length(srcData));
  for i := 0 to High(data) do
    data[i] := (srcData[i] * (1 shl (16 - bitRange))) / High(SmallInt);
  dct := TEncoder.ComputeDCT(Length(data), data);
end;

procedure TChunk.InitMix;
begin
  mixList.Clear;
end;

procedure TChunk.AddToMix(chk: TChunk);
begin
  mixList.Add(chk);
end;

procedure TChunk.FinalizeMix;
var
  i, k: Integer;
  acc: Double;
begin
  if not Assigned(mixList) or (mixList.Count = 0) then
    Exit;

  SetLength(dstData, band.frame.encoder.chunkSize);

  for k := 0 to band.frame.encoder.chunkSize - 1 do
  begin
    acc := 0.0;

    for i := 0 to mixList.Count - 1 do
      acc += mixList[i].srcData[k] * (1 shl (16 - mixList[i].bitRange));

    dstData[k] := TEncoder.makeOutputSample(acc / mixList.Count, band.frame.encoder.StoredBitDepth, band.frame.encoder.StoredBitDepth);
  end;
end;

procedure TChunk.ComputeBitRange;
var
  i, hiSmp: Integer;
begin
  hiSmp := 0;
  for i := 0 to band.frame.encoder.chunkSize - 1 do
    hiSmp := max(hiSmp, abs(TEncoder.make16BitSample(srcData[i])));
  bitRange := EnsureRange(1 + ceil(log2(hiSmp + 1.0)), 24 - band.frame.encoder.OutputBitDepth, 16);
  dstBitShift := bitRange - band.frame.encoder.StoredBitDepth;
end;

procedure TChunk.MakeSrcData(origData: PDouble);
var
  j, k, pos, n: Integer;
  acc: Double;
begin
  for j := 0 to band.frame.encoder.chunkSize - 1 do
  begin
    pos := j * band.globalData^.underSample;

    acc := 0.0;
    n := 0;
    for k := 0 to band.globalData^.underSample - 1 do
    begin
      if pos + k >= band.frame.SampleCount then
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
  SetLength(dstData, band.frame.encoder.chunkSize);
  for i := 0 to band.frame.encoder.chunkSize - 1 do
    dstData[i] := TEncoder.makeOutputSample(srcData[i], band.frame.encoder.StoredBitDepth, dstBitShift);
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
var
  i, j, k, l: Integer;
  a, b : Integer;
  cl: TChunkList;
begin
  for i := 0 to FrameCount - 1 do
  begin
    for j := 0 to BandCount - 1 do
    begin
      Assert(frames[i].bands[j].chunks.Count mod (BlockSampleCount div ChunkSize) = 0);
      AStream.WriteWord(frames[i].bands[j].chunks.Count div (BlockSampleCount div ChunkSize));

      cl := frames[i].bands[j].reducedChunks;
      if cl.Count = 0 then
        cl := frames[i].bands[j].chunks;

      for k := 0 to cl.Count - 1 do
        for l := 0 to ChunkSize - 1 do
        begin
          if l < Length(cl[k].dstData) then
            AStream.WriteByte(cl[k].dstData[l])
          else
            AStream.WriteByte($80);
        end;
    end;

    for j := 0 to BandCount - 1 do
    begin
      cl := frames[i].bands[j].chunks;

      for k := 0 to cl.Count - 1 do
      begin
        if (k and 7 = 0) and (OutputBitDepth > 8) then
        begin
          b := 0;
          for l := k to k + 7 do
          begin
            b := b shl 1;
            if l < cl.Count then
              b := b or (cl[l].dstBitShift shr 3);
          end;
          AStream.WriteByte(b);
        end;

        a := cl[k].reducedChunk.index;
        AStream.WriteByte(a and $ff);
      end;
    end;
  end;

  AStream.WriteWord(0); // Termination
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
  hc, ratioP, ratioL: Double;
  bnd: TBandGlobalData;
begin
  for i := 0 to BandCount - 1 do
  begin
    FillChar(bnd, SizeOf(bnd), 0);

    // determing low and high bandpass frequencies

    hc := min(HighCut, SampleRate / 2);
    ratioP := round((BandBoundsOffset - log2(hc)) / BandCount);
    ratioL := 0.5;//hc / SampleRate;

    if i = 0 then
      bnd.fcl := LowCut / SampleRate
    else
      bnd.fcl := power(2.0, (BandCount - i) * ratioP) * ratioL;

    if i = BandCount - 1 then
      bnd.fch := hc / SampleRate
    else
      bnd.fch := power(2.0, (BandCount - 1 - i) * ratioP) * ratioL;

    // undersample if the band high freq is a lot lower than nyquist

    bnd.underSample := Max(1, round(0.25 / bnd.fch));

    bandData[i] := bnd;
  end;
end;

procedure TEncoder.MakeBandSrcData(AIndex: Integer);
var
  j: Integer;
  bnd: TBandGlobalData;
  cic: TCICFilter;
  cicSmp: Double;
begin
  bnd := bandData[AIndex];

  SetLength(bnd.filteredData, SampleCount);
  for j := 0 to SampleCount - 1 do
    bnd.filteredData[j] := makeFloatSample(srcData[j]);

  // band pass the samples
  bnd.filteredData := DoBPFilter(bnd.fcl, bnd.fch, BandTransFactor, bnd.filteredData);
  SetLength(bnd.filteredData, SampleCount);

  if (AIndex < BandCount - 1) or TrebleBoost then
  begin
    // compensate for decoder altering the pass band (low pass)
    cic := TCICFilter.Create(bnd.fch, BandDenoiseStages, True);
    try
      for j := -cic.Ratio + 1 to High(bnd.filteredData) do
      begin
        cicSmp := cic.ProcessSample(bnd.filteredData[Min(High(bnd.filteredData), j + cic.Ratio - 1)]);
        if j >= 0 then
          bnd.filteredData[j] += 2.4 * cicSmp;
      end;
    finally
      cic.Free;
    end;
  end;

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
  BlockSampleCount := Max(trunc(8 / log2(OutputBitDepth - 7)) * ChunkSize, BlockSampleCount);

  // ensure srcData ends on a full block
  psc := SampleCount;
  SampleCount := ((SampleCount - 1) div BlockSampleCount + 1) * BlockSampleCount;
  SetLength(srcData, SampleCount);
  for i :=psc to SampleCount - 1 do
    srcData[i] := 0;

  ProjectedByteSize := ceil((SampleCount / SampleRate) * BitRate * (1024 / 8));

  fixedCost := 0;
  frameCost := 0;
  for i := 0 to BandCount - 1 do
  begin
    fixedCost += (SampleCount * (round(log2(ChunksPerBand)) + round(log2(OutputBitDepth - 7)))) div (8 * ChunkSize * bandData[i].underSample) + SizeOf(Word);
    frameCost += ChunksPerBand * ChunkSize;
  end;

  FrameCount := (ProjectedByteSize - fixedCost - 1) div frameCost + 1;

  ProjectedByteSize := fixedCost + FrameCount * frameCost;

  writeln('SampleRate = ', SampleRate, ifthen(SampleRate <> MDSampleRate, ' /!\ Won''t play properly on MegaDrive'));
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
    avgPower += sqr(makeFloatSample(srcData[i]));
  avgPower /= BandCount * SampleCount;

  totalPower := 0.0;
  for i := 0 to SampleCount - 1 do
    totalPower += avgPower - avgPower * VariableFrameSizeRatio + sqr(VariableFrameSizeRatio * makeFloatSample(srcData[i]));

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
    curPower += avgPower - avgPower * VariableFrameSizeRatio + sqr(VariableFrameSizeRatio * makeFloatSample(srcData[i]));

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

  for i := 0 to BandCount - 1 do
     WriteLn('Band #', i, ' (', round(bandData[i].fcl * SampleRate), ' Hz .. ', round(bandData[i].fch * SampleRate), ' Hz); ', bandData[i].underSample);
end;

procedure TEncoder.MakeFrames;

  procedure DoFrame(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    frm: TFrame;
  begin
    frm := frames[AIndex];

    frm.MakeBands;
  end;
begin
  ProcThreadPool.DoParallelLocalProc(@DoFrame, 0, FrameCount - 1, nil);
  WriteLn;
end;

function TEncoder.DoFilter(const samples, coeffs: TDoubleDynArray): TDoubleDynArray;
var
  i: Integer;
begin
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
    coeffs := DoFilterCoeffs(fcl, transFactor * (exp(fcl) - 1.0), True, True);
    Result := DoFilter(Result, coeffs);
  end;

  if fch < 0.5 then
  begin
    coeffs := DoFilterCoeffs(fch, transFactor * (exp(fch) - 1.0), False, True);
    Result := DoFilter(Result, coeffs);
  end;
end;

constructor TEncoder.Create(InFN, OutFN: String);
begin
  inputFN := InFN;
  outputFN := OutFN;

  BitRate := 128;
  Precision := 7;
  LowCut := 32.0;
  HighCut := 18000.0;
  OutputBitDepth := 9;
  ChunkSize := 4;
  AlternateReduce := False;
  TrebleBoost := False;
  VariableFrameSizeRatio := 0.5;

  StoredBitDepth := 8;
  ChunksPerBand := 256;
  BandBoundsOffset := 8;
  BandDenoiseStages := 2;
  BandTransFactor := 0.1;

  frames := TFrameList.Create;
end;

destructor TEncoder.Destroy;
begin
  frames.Free;

  inherited Destroy;
end;

procedure TEncoder.MakeDstData;
var
  i, j, k, pos, poso: Integer;
  smp: Double;
  offset: array[0 .. BandCount - 1] of Integer;
  bnd: TBandGlobalData;
  lp: array[0 .. BandCount - 1] of TCICFilter;
begin
  WriteLn('MakeDstData');

  SetLength(dstData, SampleCount);
  FillWord(dstData[0], Length(dstData), 0);

  for i := 0 to BandCount - 1 do
  begin
    bnd := bandData[i];

    SetLength(bnd.dstData, Length(dstData));
    FillWord(bnd.dstData[0], Length(dstData), 0);

    offset[i] := 0;
    lp[i] := Nil;

    if i < BandCount - 1 then
    begin
      lp[i] := TCICFilter.Create(bnd.fch, BandDenoiseStages, False);
      offset[i] += -lp[i].Ratio + 1;
    end;

    bandData[i] := bnd;
  end;

  pos := 0;
  try
    for k := 0 to frames.Count - 1 do
      for i := 0 to frames[k].SampleCount - 1 do
      begin
        for j := 0 to BandCount - 1 do
        begin
          smp := frames[k].bands[j].dstData[i];

          if Assigned(lp[j]) then
            smp := lp[j].ProcessSample(smp);

          poso := pos + offset[j];
          if InRange(poso, 0, High(dstData)) then
          begin
            bandData[j].dstData[poso] := make16BitSample(smp);
            dstData[poso] := EnsureRange(dstData[poso] + bandData[j].dstData[poso], Low(SmallInt), High(SmallInt));
          end;
        end;

        Inc(pos);
      end;
  finally
    for i := 0 to BandCount - 1 do
      lp[i].Free;
  end;
end;

function TEncoder.DoFilterCoeffs(fc, transFactor: Double; HighPass, Windowed: Boolean): TDoubleDynArray;
var
  sinc, win, sum: Double;
  i, N: Integer;
begin
  N := ceil(4.6 / transFactor);
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
  smp16: SmallInt;
begin
  smp16 := round(smp * (1 shl (15 - bitShift)));
  Result := EnsureRange(smp16 + (1 shl (OutBitDepth - 1)), 0, (1 shl OutBitDepth) - 1);
end;

class function TEncoder.makeFloatSample(smp: SmallInt): Double;
begin
  Result := smp / -Low(SmallInt);
end;

class function TEncoder.makeFloatSample(smp: Byte; OutBitDepth, bitShift: Integer): Double;
var
  smp16: SmallInt;
begin
  smp16 := SmallInt(smp) - (1 shl (OutBitDepth - 1));
  Result := smp16 / (1 shl (15 - bitShift));
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

class function TEncoder.CompareDCT(firstCoeff, lastCoeff: Integer; compress: Boolean; const dctA, dctB: TDoubleDynArray): Double;
var
  i: Integer;
begin
  Result := 0.0;
  if compress then
  begin
    for i := firstCoeff to lastCoeff do
      Result += sqr(CompressDCT(dctA[i]) - CompressDCT(dctB[i]));
  end
  else
  begin
    for i := firstCoeff to lastCoeff do
      Result += sqr(dctA[i] - dctB[i]);
  end;
  Result := sqrt(Result);
end;

class function TEncoder.CompressDCT(coeff: Double): Double;
begin
  Result := Sign(coeff) * power(Abs(coeff), 0.707);
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

var
  enc: TEncoder;
  i: Integer;
begin
  try
    FormatSettings.DecimalSeparator := '.';

{$ifdef DEBUG}
    ProcThreadPool.MaxThreadCount := 1;
{$endif}

    if ParamCount < 2 then
    begin
      WriteLn('Usage: ', ExtractFileName(ParamStr(0)) + ' <source file> <dest file> [options]');
      Writeln('Main options:');
      WriteLn(#9'-br'#9'encoder bit rate in kilobits/second; example: "-br128"');
      WriteLn(#9'-lc'#9'bass cutoff frequency');
      WriteLn(#9'-hc'#9'treble cutoff frequency');
      WriteLn(#9'-tb'#9'apply treble boost');
      WriteLn(#9'-vfr'#9'RMS power based variable frame size ratio (0.0-1.0)');
      WriteLn(#9'-v'#9'verbose mode');
      Writeln('Development options:');
      WriteLn(#9'-cs'#9'chunk size');
      WriteLn(#9'-ar'#9'use alternate clustering reduce method');
      WriteLn(#9'-obd'#9'output bit depth (8-16)');
      WriteLn(#9'-pr'#9'K-means precision; 0: "lossless" mode');

      WriteLn;
      Writeln('(source file must be 16bit mono WAV, preferably 26390Hz for MegaDrive playback)');
      WriteLn;
      Exit;
    end;

    enc := TEncoder.Create(ParamStr(1), ParamStr(2));
    try
      enc.BitRate := round(ParamValue('-br', enc.BitRate));
      enc.Precision := round(ParamValue('-pr', enc.Precision));
      enc.LowCut :=  ParamValue('-lc', enc.LowCut);
      enc.HighCut :=  ParamValue('-hc', enc.HighCut);
      enc.TrebleBoost := HasParam('-tb');
      enc.VariableFrameSizeRatio :=  ParamValue('-vfr', enc.VariableFrameSizeRatio);
      enc.OutputBitDepth :=  round(ParamValue('-obd', enc.OutputBitDepth));
      enc.ChunkSize :=  round(ParamValue('-cs', enc.ChunkSize));
      enc.Verbose := HasParam('-v');
      enc.AlternateReduce := HasParam('-ar');

      WriteLn('BitRate = ', FloatToStr(enc.BitRate));
      WriteLn('LowCut = ', FloatToStr(enc.LowCut));
      WriteLn('HighCut = ', FloatToStr(enc.HighCut));
      WriteLn('TrebleBoost = ', BoolToStr(enc.TrebleBoost, True));
      WriteLn('VariableFrameSizeRatio = ', FloatToStr(enc.VariableFrameSizeRatio));
      if enc.Verbose then
      begin
        WriteLn('ChunkSize = ', enc.ChunkSize);
        WriteLn('AlternateReduce = ', BoolToStr(enc.AlternateReduce, True));
        WriteLn('OutputBitDepth = ', enc.OutputBitDepth);
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

      enc.ComputeEAQUAL(enc.SampleCount, True, True, enc.srcData, enc.dstData);

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

