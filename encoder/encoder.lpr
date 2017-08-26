program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, yakmo, ap, fft, conv, anysort, minlbfgs, kmeans, correlation;

const
  BandCount = 2;
  InputSNRDb = -90.3; // 16bit

type
  TEncoder = class;
  TFrame = class;
  TBand = class;
  TChunk = class;

  TBandGlobalData = record
    fcl, fch: Double;
    underSample: Integer;
    chunkSize: Integer;
    chunkCount: Integer;
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

  { TSimpleFilter }

  TSimpleFilter = class
  private
    prevSmp: TDoubleDynArray;
  public
    Count: Integer;
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

    procedure ComputeDCT(FromDstData: Boolean);
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
    rmsPower: Double;
    desiredChunkCount: Integer;
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
    procedure MakeDstData;
  end;

  { TFrame }

  TFrame = class
  public
    encoder: TEncoder;

    index: Integer;
    SampleCount: Integer;
    FrameSize: Integer;
    BandRMSWeighting: Boolean;
    BandBWeighting: Boolean; // otherwise A-weighting
    BandWeightingApplyPower: Double;

    bands: array[0..BandCount - 1] of TBand;

    constructor Create(enc: TEncoder; idx: Integer; startSample, endSample: Integer);
    destructor Destroy; override;

    procedure FindBandDesiredChunkCount;
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
    BandDealiasSecondOrder: Boolean; // otherwise first order
    StoredBitDepth: Integer; // max 8Bits
    OutputBitDepth: Integer; // 8 to 16Bits
    MinChunkSize: Integer;
    MaxFrameSize: Integer;
    MaxChunksPerBand: Integer;
    AlternateReduce: Boolean;
    BandBoundsOffset: Integer;
    BandDenoiseStages: Integer;

    SampleRate: Integer;
    SampleCount: Integer;
    projectedByteSize, frameCount, frameSampleCount: Integer;
    verbose: Boolean;

    BWSearch, CRSearch: Boolean;

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

    procedure SearchBestBandWeighting;

    function DoFilterCoeffs(fc, transFactor: Double; HighPass, Windowed: Boolean): TDoubleDynArray;
    function DoFilter(const samples, coeffs: TDoubleDynArray): TDoubleDynArray;
    function DoBPFilter(fcl, fch, transFactor: Double; const samples: TDoubleDynArray): TDoubleDynArray;

    function ComputeEAQUAL(chunkSz: Integer; UseDIX, Verbz: Boolean; const smpRef, smpTst: TSmallIntDynArray): Double;
  end;


function IsDebuggerPresent () : LongBool stdcall; external 'kernel32.dll';

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

constructor TSimpleFilter.Create(AFc: Double; AStages: Integer; AHighPass: Boolean);
begin
  Stages := AStages;
  HighPass := AHighPass;
  Count := trunc(sqrt(0.196202 + sqr(AFc)) / AFc);
  SetLength(prevSmp, Stages);
end;

function TSimpleFilter.ProcessSample(Smp: Double): Double;
var
  i: Integer;
  v: Double;
begin
  Result := Smp;
  for i := 0 to Stages - 1 do
  begin
    v := prevSmp[i];
    prevSmp[i] := Result;
    Result := (Count * v + Result) / (Count + 1);
  end;

  if HighPass then
    Result := Smp - Result;
end;

constructor TFrame.Create(enc: TEncoder; idx: Integer; startSample, endSample: Integer);
var
  i: Integer;
begin
  encoder := enc;
  index := idx;

  BandRMSWeighting := True;
  BandBWeighting := False;
  BandWeightingApplyPower := 0.0;

  SampleCount := encoder.frameSampleCount;

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


procedure TFrame.FindBandDesiredChunkCount;
var
  i, j, allSz, dbsOh: Integer;
  sz, dbh: Double;
  wgt, fsq: Double;
  full: Boolean;
  bnd: TBand;
  s: String;
  dcc: array[0..BandCount - 1] of Integer;
begin
  dbh := -MaxDouble;
  for i := 0 to BandCount - 1 do
    dbh := max(dbh, bands[i].rmsPower);
  dbh := max(InputSNRDb + 1.0, dbh);

  for i := 0 to BandCount - 1 do
  begin
    bnd := bands[i];

    fsq := bnd.globalData^.fcl * bnd.globalData^.fch * encoder.SampleRate * encoder.SampleRate;

    if BandBWeighting then
      // B-weighting
      wgt := sqr(12194.0) * sqrt(fsq) * fsq / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)) * sqrt(fsq + sqr(158.5)))
    else
      // A-weighting
      wgt := sqr(12194.0) * sqr(fsq) / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)) * sqrt((fsq + sqr(107.7)) * (fsq + sqr(737.9))));

    bnd.weight := power(wgt, BandWeightingApplyPower);
    if BandRMSWeighting then
      bnd.weight *= max(1.0, (bnd.rmsPower - InputSNRDb)) / (dbh - InputSNRDb);
  end;

  dbsOh := 0;
  for i := 0 to BandCount - 1 do
  begin
    bnd := bands[i];

    Inc(dbsOh);
    for j := 1 to bnd.chunks.Count - 1 do
      if bnd.chunks[j].dstBitShift <> bnd.chunks[j - 1].dstBitShift then
        Inc(dbsOh);
  end;

  allSz := 0;
  sz := 1;
  repeat
    full := True;
    allSz := dbsOh * SizeOf(Byte) + SizeOf(Word);
    for i := 0 to BandCount - 1 do
    begin
      bnd := bands[i];

      dcc[i] := round(sz * bnd.weight);
      dcc[i] := EnsureRange(dcc[i], 1, min(encoder.MaxChunksPerBand, bnd.globalData^.chunkCount));

      full := full and not (dcc[i] < min(encoder.MaxChunksPerBand, bnd.globalData^.chunkCount));

      allSz += dcc[i] * bnd.globalData^.chunkSize + bnd.globalData^.chunkCount * SizeOf(Byte) + SizeOf(Word);
    end;

    if allSz <= encoder.MaxFrameSize then
    begin
      FrameSize := allSz;
      for i := 0 to BandCount - 1 do
        bands[i].desiredChunkCount := dcc[i];
    end;

    sz += 0.1;
  until full;

  if encoder.verbose and not (encoder.BWSearch or encoder.CRSearch) then
  begin
    s := '#' + IntToStr(index) + ', ' + IntToStr(FrameSize) + #9;
    for i := 0 to BandCount - 1 do
      s := s + IntToStr(bands[i].desiredChunkCount) + ', ' + FormatFloat('0.0', bands[i].rmsPower) + #9;
    WriteLn(s);;
  end;
end;

procedure TFrame.MakeBands;
var
  i: Integer;
begin
  for i := 0 to BandCount - 1 do
    bands[i].MakeChunks;

  FindBandDesiredChunkCount;

  for i := 0 to BandCount - 1 do
  begin
    if not encoder.CRSearch then
      bands[i].KMeansReduce;

    bands[i].MakeDstData;
  end;

  if not encoder.verbose or encoder.BWSearch then
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
var i: Integer;
begin
  frame := frm;
  index := idx;
  globalData := @frame.encoder.bandData[index];
  AlternateReduceMethod := frame.encoder.AlternateReduce;

  srcData := @globalData^.filteredData[startSample];

  rmsPower := 0.0;
  for i := startSample to endSample do
    rmsPower += sqr(globalData^.filteredData[i]);
  rmsPower := sqrt(rmsPower / (endSample - startSample + 1)) - 1.0 / Low(SmallInt);
  rmsPower := 20.0 * log10(rmsPower);

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
  chunks.Capacity := globalData^.chunkCount;
  for i := 0 to globalData^.chunkCount - 1 do
  begin
    chunk := TChunk.Create(Self, i);
    chunk.ComputeBitRange;
    chunk.MakeDstData;
    if not frame.encoder.CRSearch then
      chunk.ComputeDCT(False);
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
  if (prec = 0) or (desiredChunkCount = globalData^.chunkCount) then Exit;

  //WriteLn('KMeansReduce #', index, ' ', globalData^.desiredChunkCount);

  SetLength(XY, chunks.Count, globalData^.chunkSize);

  for i := 0 to chunks.Count - 1 do
  begin
    //// add first and last sample to features to allow continuity between chunks
    //XY[i, 0] := chunks[i].srcData[0];
    //XY[i, 1] := chunks[i].srcData[globalData^.chunkSize - 1];
    //
    //// add approximate derivatives of first and last sample to features to allow continuity between chunks
    //XY[i, 2] := -1.5 * chunks[i].srcData[0] +
    //        2.0 * chunks[i].srcData[1] +
    //        -0.5 * chunks[i].srcData[2];
    //XY[i, 3] := 1.5 * chunks[i].srcData[globalData^.chunkSize - 1] +
    //        -2.0 * chunks[i].srcData[globalData^.chunkSize - 2] +
    //        0.5 * chunks[i].srcData[globalData^.chunkSize - 3];
    //
    for j := 0 to globalData^.chunkSize - 1 do
      XY[i, j] := chunks[i].dct[j];
  end;

  DoExternalKMeans(XY, desiredChunkCount, 1, prec, AlternateReduceMethod, False, XYC);

  reducedChunks.Clear;
  reducedChunks.Capacity := desiredChunkCount;
  for i := 0 to desiredChunkCount - 1 do
  begin
    chunk := TChunk.Create(Self, i);

    reducedChunks.Add(chunk);
    DoXYC(i);
  end;
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

    Assert((frame.encoder.Precision = 0) or (globalData^.chunkCount = desiredChunkCount) or (chunk.reducedChunk <> chunk));

    for j := 0 to globalData^.chunkSize - 1 do
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

  SetLength(srcData, band.globalData^.chunkSize);

  origSrcData := @band.srcData[idx * band.globalData^.chunkSize * band.globalData^.underSample];

  MakeSrcData(origSrcData);

  reducedChunk := Self;

  mixList := TChunkList.Create(False);
end;

destructor TChunk.Destroy;
begin
  mixList.Free;

  inherited Destroy;
end;

procedure TChunk.ComputeDCT(FromDstData: Boolean);
var
  i: Integer;
  data: TDoubleDynArray;
begin
  SetLength(data, band.globalData^.chunkSize);
  if FromDstData then
  begin
    for i := 0 to High(data) do
      data[i] := (Integer(dstData[i]) + Low(ShortInt)) / High(ShortInt);
  end
  else
  begin
    for i := 0 to High(data) do
      data[i] := (srcData[i] * (1 shl (16 - bitRange))) / High(SmallInt);
  end;

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

  SetLength(dstData, band.globalData^.chunkSize);

  for k := 0 to band.globalData^.chunkSize - 1 do
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
  for i := 0 to band.globalData^.chunkSize - 1 do
    hiSmp := max(hiSmp, abs(TEncoder.make16BitSample(srcData[i])));
  bitRange := EnsureRange(1 + ceil(log2(hiSmp + 1.0)), 24 - band.frame.encoder.OutputBitDepth, 16);
  dstBitShift := bitRange - band.frame.encoder.StoredBitDepth;
end;

procedure TChunk.MakeSrcData(origData: PDouble);
var
  j, k, pos, n: Integer;
  acc: Double;
begin
  for j := 0 to band.globalData^.chunkSize - 1 do
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
  SetLength(dstData, band.globalData^.chunkSize);
  for i := 0 to band.globalData^.chunkSize - 1 do
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
  writeln('SampleRate = ', SampleRate);
  WriteLn('SampleCount = ', SampleCount);
end;

procedure TEncoder.SaveWAV;
var
  fs: TFileStream;
  wavFN: String;
begin
  wavFN := ChangeFileExt(outputFN, '.wav');

  if not CRSearch then
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
  i, j, k: Integer;
  ms: TMemoryStream;
  b, pb : Integer;
  cl: TChunkList;
begin
  AStream.WriteWord(bandData[0].chunkCount);

  ms := TMemoryStream.Create;
  try
    for i := 0 to frameCount - 1 do
    begin
      for j := 0 to BandCount - 1 do
      begin
        cl := frames[i].bands[j].reducedChunks;
        if cl.Count = 0 then
          cl := frames[i].bands[j].chunks;
        ms.WriteWord(cl.Count);

        for k := 0 to cl.Count - 1 do
        begin
          ms.WriteBuffer(cl[k].dstData[0], bandData[j].chunkSize);
        end;
      end;

      for j := 0 to BandCount - 1 do
      begin
        pb := -1;
        for k := 0 to bandData[j].chunkCount - 1 do
        begin
          b := frames[i].bands[j].chunks[k].dstBitShift;
          if b <> pb then
          begin
            ms.WriteByte(255 - StoredBitDepth + b);
            pb := b;
          end;

          ms.WriteByte(frames[i].bands[j].chunks[k].reducedChunk.index);
        end;
      end;

      ms.Position := 0;
      AStream.WriteWord(ms.Size);
      AStream.CopyFrom(ms, ms.Size);
      ms.Clear;
    end;
  finally
    ms.Free;
  end;
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

    bnd.chunkSize := minChunkSize;

    // undersample if the band high freq is a lot lower than nyquist

    bnd.underSample := Max(1, round(0.125 / bnd.fch));

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
  SetLength(bnd.filteredData, frameSampleCount * frameCount);

  for j := SampleCount to frameSampleCount * frameCount - 1 do
    bnd.filteredData[j] := 0;

  if AIndex < BandCount - 1 then
  begin
    // compensate for decoder altering the pass band (low pass)
    cic := TCICFilter.Create(bnd.fch, BandDenoiseStages, False);
    try
      for j := -cic.Ratio + 1 to High(bnd.filteredData) do
      begin
        cicSmp := cic.ProcessSample(bnd.filteredData[Min(High(bnd.filteredData), j + cic.Ratio - 1)]);
        if j >= 0 then
          bnd.filteredData[j] += cic.Stages * (bnd.filteredData[j]- cicSmp);
      end;
    finally
      cic.Free;
    end;
  end;

  if AIndex > 0 then
  begin
    // compensate for decoder altering the pass band (high pass)
    cic := TCICFilter.Create(bnd.fcl, BandDenoiseStages, False);
    try
      for j := -cic.Ratio + 1 to High(bnd.filteredData) do
      begin
        cicSmp := cic.ProcessSample(bnd.filteredData[Min(High(bnd.filteredData), j + cic.Ratio - 1)]);
        if j >= 0 then
          bnd.filteredData[j] += cicSmp;
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
  i, blockSampleCount, blockChunkCount, curStart, curEnd: Integer;
  frm: TFrame;
begin
  MaxChunksPerBand := 255 - (OutputBitDepth - StoredBitDepth);

  MakeBandGlobalData;

  blockSampleCount := 0;
  for i := 0 to BandCount - 1 do
    if bandData[i].underSample > blockSampleCount then
      blockSampleCount := bandData[i].underSample * bandData[i].chunkSize;

  blockChunkCount := 0;
  for i := 0 to BandCount - 1 do
    blockChunkCount += blockSampleCount div (bandData[i].underSample * bandData[i].chunkSize);

  // pass 1
  projectedByteSize := ceil((SampleCount / SampleRate) * BitRate * (1024 / 8));
  frameCount :=  projectedByteSize div MaxFrameSize;
  frameSampleCount := (SampleCount - 1) div frameCount + 1;
  frameSampleCount := ((frameSampleCount - 1) div blockSampleCount + 1) * blockSampleCount;

  // pass 2
  projectedByteSize := MaxFrameSize * frameCount;

  for i := 0 to BandCount - 1 do
    bandData[i].chunkCount := (frameSampleCount - 1) div (bandData[i].chunkSize * bandData[i].underSample) + 1;

  if verbose and not (BWSearch or CRSearch) then
  begin
    writeln('MaxFrameSize = ', MaxFrameSize);
    writeln('FrameCount = ', frameCount);
    writeln('FrameSampleCount = ', frameSampleCount);
    writeln('ProjectedByteSize = ', projectedByteSize);
  end;

  if not (BWSearch or CRSearch) then
    for i := 0 to BandCount - 1 do
       WriteLn('Band #', i, ' (', round(bandData[i].fcl * SampleRate), ' Hz .. ', round(bandData[i].fch * SampleRate), ' Hz); ', bandData[i].chunkSize, ' * ', bandData[i].chunkCount, '; ', bandData[i].underSample);

  ProcThreadPool.DoParallelLocalProc(@DoBand, 0, BandCount - 1, nil);

  frames.Clear;
  for i := 0 to frameCount - 1 do
  begin
    curStart := i * frameSampleCount;
    curEnd := min(SampleCount, (i + 1) * frameSampleCount) - 1;

    frm := TFrame.Create(Self, i, curStart, curEnd);
    frames.Add(frm);
  end;
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
  ProcThreadPool.DoParallelLocalProc(@DoFrame, 0, frameCount - 1, nil);
  if not CRSearch then
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
  BandTransFactor := 0.1;
  LowCut := 32.0;
  HighCut := 18000.0;
  BandDealiasSecondOrder := True;
  StoredBitDepth := 8;
  OutputBitDepth := 13;
  MinChunkSize := 8;
  MaxFrameSize:= 16 * 1024;
  AlternateReduce := False;

  BandBoundsOffset := 8;
  BandDenoiseStages := 2;

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
  lp, hp: array[0 .. BandCount - 1] of TCICFilter;
begin
  if not (BWSearch or CRSearch) then
    WriteLn('MakeDstData');

  SetLength(dstData, frameSampleCount * frames.Count);
  FillWord(dstData[0], Length(dstData), 0);

  for i := 0 to BandCount - 1 do
  begin
    bnd := bandData[i];

    SetLength(bnd.dstData, Length(dstData));
    FillWord(bnd.dstData[0], Length(dstData), 0);

    offset[i] := 0;
    lp[i] := Nil;
    hp[i] := Nil;

    if i < BandCount - 1 then
    begin
      lp[i] := TCICFilter.Create(bnd.fch, BandDenoiseStages, False);
      offset[i] += -lp[i].Ratio + 1;
    end;

    if i > 0 then
    begin
      hp[i] := TCICFilter.Create(bnd.fcl, BandDenoiseStages, True);
      offset[i] += -hp[i].Ratio + 1;
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

          if Assigned(hp[j]) then
            smp := hp[j].ProcessSample(smp);

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
    begin
      lp[i].Free;
      hp[i].Free;
    end;
  end;
end;

procedure TEncoder.SearchBestBandWeighting;
type
  TSBW = record
    RMSWeighting: Boolean;
    BWeighting: Boolean;
    ApplyPower: Integer;
    AltMethod: Integer;
  end;

const
  sbw: array[0..3] of TSBW =
  (
    // per band
    (RMSWeighting: True;  BWeighting: True;  ApplyPower: 0;  AltMethod: 0),
    (RMSWeighting: True;  BWeighting: True;  ApplyPower: 0;  AltMethod: 1),
    (RMSWeighting: True;  BWeighting: True;  ApplyPower: 0;  AltMethod: 2),
    // reference
    (RMSWeighting: True;  BWeighting: True;  ApplyPower: 0;  AltMethod: -1)
    // combinations
    //(RMSWeighting: True;  BWeighting: True;  ApplyPower: 1;  AltMethod: -1),
    //(RMSWeighting: True;  BWeighting: True;  ApplyPower: -1; AltMethod: -1),
    //(RMSWeighting: True;  BWeighting: False; ApplyPower: 1;  AltMethod: -1),
    //(RMSWeighting: True;  BWeighting: False; ApplyPower: -1; AltMethod: -1),
    //(RMSWeighting: True;  BWeighting: True;  ApplyPower: 2;  AltMethod: -1),
    //(RMSWeighting: True;  BWeighting: True;  ApplyPower: -2; AltMethod: -1),
    //(RMSWeighting: False; BWeighting: True;  ApplyPower: 0;  AltMethod: -1),
    //(RMSWeighting: False; BWeighting: True;  ApplyPower: 1;  AltMethod: -1),
    //(RMSWeighting: False; BWeighting: True;  ApplyPower: 2;  AltMethod: -1),
    //(RMSWeighting: False; BWeighting: False; ApplyPower: 1;  AltMethod: -1)
  );

var
  results: array[0..High(sbw)] of TDoubleDynArray;

  procedure DoEAQ(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    ref, tst: TSmallIntDynArray;
  begin
    ref := Copy(srcData, AIndex * frameSampleCount, frameSampleCount);
    SetLength(ref, frameSampleCount);
    tst := Copy(dstData, AIndex * frameSampleCount, frameSampleCount);
    results[PtrInt(AData), AIndex] := ComputeEAQUAL(frameSampleCount, True, False, ref, tst);
  end;

var
  i, j, k: Integer;
  best: Double;
begin
  BWSearch := True;
  try
    for i := 0 to High(sbw) do
    begin
      PrepareFrames;

      for j := 0 to frameCount - 1 do
      begin
        frames[j].BandRMSWeighting := sbw[i].RMSWeighting;
        frames[j].BandBWeighting := sbw[i].BWeighting;
        frames[j].BandWeightingApplyPower := sbw[i].ApplyPower;
        for k := 0 to BandCount - 1 do
          frames[j].bands[k].AlternateReduceMethod := sbw[i].AltMethod = k;
      end;

      MakeFrames;
      MakeDstData;

      SetLength(results[i], frameCount);

      ProcThreadPool.DoParallelLocalProc(@DoEAQ, 0, frameCount - 1, Pointer(i));
    end;
  finally
    BWSearch := False;
  end;

  PrepareFrames;

  for j := 0 to frameCount - 1 do
  begin
    for i := 0 to BandCount - 1 do
      frames[j].bands[i].AlternateReduceMethod := results[i, j] > results[BandCount, j];

    best := -MaxDouble;
    for i := BandCount to High(sbw) do
      if results[i, j] > best then
      begin
        best := results[i, j];
        frames[j].BandRMSWeighting := sbw[i].RMSWeighting;
        frames[j].BandBWeighting := sbw[i].BWeighting;
        frames[j].BandWeightingApplyPower := sbw[i].ApplyPower;
      end;

    Write('#', j, #9, frames[j].BandRMSWeighting, #9, frames[j].BandBWeighting, #9, FloatToStr(frames[j].BandWeightingApplyPower));
    for i := 0 to BandCount - 1 do
      Write(#9, frames[j].bands[i].AlternateReduceMethod);
    WriteLn;
  end;

  MakeFrames;
  MakeDstData;
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
      WriteLn(#9'-br'#9'encoder bit rate in kilobits/second; example: "-q128"');
      WriteLn(#9'-pr'#9'K-means precision; 0: "lossless" mode');
      WriteLn(#9'-lc'#9'bass cutoff frequency');
      WriteLn(#9'-hc'#9'treble cutoff frequency');
      WriteLn(#9'-btf'#9'band transition factor (0.0001-1)');
      WriteLn(#9'-bd1'#9'use first order dealias filter (otherwise second order)');
      WriteLn(#9'-obd'#9'output bit depth (1-8)');
      WriteLn(#9'-mcs'#9'minimum chunk size');
      WriteLn(#9'-sbw'#9'per frame search for best band weighting');
      WriteLn(#9'-v'#9'verbose K-means');
      WriteLn(#9'-al'#9'use alternate clusering reduce method');
      WriteLn;
      Writeln('(source file must be 16bit mono WAV)');
      WriteLn;
      Exit;
    end;

    enc := TEncoder.Create(ParamStr(1), ParamStr(2));
    try
      enc.BitRate := round(ParamValue('-br', enc.BitRate));
      enc.Precision := round(ParamValue('-pr', enc.Precision));
      enc.LowCut :=  ParamValue('-lc', enc.LowCut);
      enc.HighCut :=  ParamValue('-hc', enc.HighCut);
      enc.BandTransFactor :=  ParamValue('-btf', enc.BandTransFactor);
      enc.BandDealiasSecondOrder :=  ParamStart('-bd1') = -1;
      enc.OutputBitDepth :=  round(ParamValue('-obd', enc.OutputBitDepth));
      enc.MinChunkSize :=  round(ParamValue('-mcs', enc.MinChunkSize));
      enc.verbose := ParamStart('-v') <> -1;
      enc.AlternateReduce := ParamStart('-ar') <> -1;

      WriteLn('BitRate = ', FloatToStr(enc.BitRate));
      WriteLn('Precision = ', enc.Precision);
      WriteLn('LowCut = ', FloatToStr(enc.LowCut));
      WriteLn('HighCut = ', FloatToStr(enc.HighCut));
      WriteLn('BandTransFactor = ', FloatToStr(enc.BandTransFactor));
      WriteLn('BandDealiasSecondOrder = ', BoolToStr(enc.BandDealiasSecondOrder, True));
      WriteLn('OutputBitDepth = ', enc.OutputBitDepth);
      WriteLn('MinChunkSize = ', enc.MinChunkSize);
      WriteLn('AlternateReduce = ', BoolToStr(enc.AlternateReduce, True));
      WriteLn;

      enc.Load;

      if ParamStart('-sbw') <> -1 then
      begin
        enc.SearchBestBandWeighting;
      end
      else
      begin
        enc.PrepareFrames;
        enc.MakeFrames;
        enc.MakeDstData;
      end;

      enc.SaveWAV;
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

