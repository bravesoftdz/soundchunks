program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, yakmo, ap, fft, conv, anysort, minlbfgs, kmeans, correlation;

const
  InputSNRDb = -90.3; // 16bit
  BandCount = 4;
  CICCorrRatio: array[Boolean{2nd order?}, 0..1] of Double = (
    (0.503, 0.816),
    (1.815, 1.815)
  );

type
  TComplex1DArrayList = specialize TFPGList<TComplex1DArray>;
  TDoubleDynArray2 = array of TDoubleDynArray;

  TEncoder = class;
  TFrame = class;
  TBand = class;

  TBandGlobalData = record
    fcl, fch: Double;
    weight: Double;

    underSample, underSampleUnMin: Integer;
    chunkSize, chunkSizeUnMin: Integer;
    chunkCount: Integer;

    filteredData: TDoubleDynArray;
    dstData: TSmallIntDynArray;
  end;

  PBandGlobalData = ^TBandGlobalData;

  { TSecondOrderCIC }

  TSecondOrderCIC = class
  private
    v, vv, v2: Double;
    R, pos: Integer;
    smps, vs: TDoubleDynArray;
  public
    Ratio: Integer;
    SecondOrder: Boolean;

    constructor Create(ARatio: Integer; ASecondOrder: Boolean);

    function ProcessSample(Smp: Double): Double;
  end;

  { TChunk }

  TChunk = class
  public
    band: TBand;
    reducedChunk, sameChunk: TChunk;

    index: Integer;
    tmpIndex: Integer;

    mixList: TComplex1DArrayList;

    srcData: TDoubleDynArray;
    fft: TComplex1DArray;
    dstData: TByteDynArray;

    constructor Create(bnd: TBand; idx: Integer);
    destructor Destroy; override;

    procedure ComputeFFT(FromDstData: Boolean);
    procedure InitMix;
    procedure AddToMix(afft: TComplex1DArray);
    procedure FinalizeMix(IntoDstData: Boolean);
    procedure MakeDstData;
  end;

  TChunkList = specialize TFPGObjectList<TChunk>;

  { TBand }

  TBand = class
  public
    frame: TFrame;

    index: Integer;
    rmsPower: Double;
    desiredChunkCount: Integer;

    bitRange: Integer;
    dstBitShift: Integer;

    srcData: PDouble;
    dstData: TDoubleDynArray;

    chunks: TChunkList;
    reducedChunks: TChunkList;
    globalData: PBandGlobalData;

    constructor Create(frm: TFrame; idx: Integer; startSample, endSample: Integer);
    destructor Destroy; override;

    procedure ComputeBitRange;

    procedure MakeChunks;

    procedure KMeansReduce;
    procedure MakeDstData;
  end;

  { TFrame }

  TFrame = class
  public
    encoder: TEncoder;

    index: Integer;
    sampleCount: Integer;

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
    BandBWeighting: Boolean; // otherwise A-weighting
    BandWeightingApplyPower: Double;
    OutputBitDepth: Integer; // max 8Bits
    MinChunkSize: Integer;
    ChunksPerFrame: Integer;
    UsePython: Boolean;

    sampleRate: Integer;
    sampleCount: Integer;
    frameByteSize, projectedByteSize, frameCount, frameSampleCount: Integer;
    verbose: Boolean;

    CRSearch: Boolean;
    CR1, CR2: Double;

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
    procedure MakeBandGlobalData2;
    procedure MakeBandSrcData(AIndex: Integer);

    procedure MakeFrames;
    procedure MakeDstData;

    procedure SearchBestCorrRatios;

    function DoFilterCoeffs(fc, transFactor: Double; HighPass: Boolean): TDoubleDynArray;
    function DoFilter(const samples, coeffs: TDoubleDynArray): TDoubleDynArray;
    function DoBPFilter(fcl, fch, transFactor: Double; const samples: TDoubleDynArray): TDoubleDynArray;

    function ComputeEAQUAL(chunkSz: Integer; UseDIX: Boolean; const smpRef, smpTst: TDoubleDynArray): Double;
    function ComputeEAQUALMulti(chunkSz: Integer; UseDIX: Boolean; const smpRef: TDoubleDynArray;
      smpTst: TDoubleDynArray2): TDoubleDynArray;
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

constructor TFrame.Create(enc: TEncoder; idx: Integer; startSample, endSample: Integer);
var
  i: Integer;
begin
  encoder := enc;
  index := idx;

  sampleCount := encoder.frameSampleCount;

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
  i, allCC: Integer;
  sz, dbh: Double;
  full: Boolean;
  bnd: TBand;
  s: String;
begin

  dbh := -MaxDouble;
  for i := 0 to BandCount - 1 do
    dbh := max(dbh, bands[i].rmsPower);
  dbh := max(InputSNRDb + 1.0, dbh);

  sz := 1;
  repeat
    full := True;
    allCC := 0;
    for i := 0 to BandCount - 1 do
    begin
      bnd := bands[i];

      bnd.desiredChunkCount := round(sz * bnd.globalData^.weight * max(1.0, (bnd.rmsPower - InputSNRDb)) / (dbh - InputSNRDb));

      bnd.desiredChunkCount := max(1, bnd.desiredChunkCount);
      bnd.desiredChunkCount := min(bnd.globalData^.chunkCount, bnd.desiredChunkCount);

      full := full and (bnd.desiredChunkCount = bnd.globalData^.chunkCount);

      allCC += bnd.desiredChunkCount * bands[i].globalData^.chunkSize div encoder.MinChunkSize;
    end;
    sz += 0.1;
  until (allCC >= encoder.ChunksPerFrame) or full;

  if encoder.verbose then
  begin
    s := '#' + IntToStr(index) + #9;
    for i := 0 to BandCount - 1 do
      s := s + IntToStr(bands[i].desiredChunkCount) + ', ' + FormatFloat('0.0', bands[i].rmsPower) + ', ' + IntToStr(bands[i].bitRange) + #9;
    WriteLn(s);;
  end;
end;

procedure TFrame.MakeBands;
var
  i: Integer;
  bnd: TBand;
begin
  FindBandDesiredChunkCount;

  for i := 0 to BandCount - 1 do
  begin
    bnd := bands[i];

    bnd.MakeChunks;

    if not encoder.CRSearch then
      bnd.KMeansReduce;

    bnd.MakeDstData;
  end;

  if not encoder.verbose then
    Write('.');
end;

constructor TSecondOrderCIC.Create(ARatio: Integer; ASecondOrder: Boolean);
begin
  SecondOrder := ASecondOrder;
  Ratio := ARatio;
  R := max(1, ARatio);
  SetLength(vs, R);
  SetLength(smps, R);
  v := 0;
  v2 := 0;
  pos := 0;
end;

function TSecondOrderCIC.ProcessSample(Smp: Double): Double;
begin
  Smp /= R;

  v := v + smp - smps[pos];
  smps[pos] := smp;

  vv := v / R;
  v2 := v2 + vv - vs[pos];
  vs[pos] := vv;

  pos := (pos + 1) mod R;

  Result := ifthen(SecondOrder, v2, v);
end;


constructor TBand.Create(frm: TFrame; idx: Integer; startSample, endSample: Integer);
var i: Integer;
begin
  frame := frm;
  index := idx;
  globalData := @frame.encoder.bandData[index];

  srcData := @globalData^.filteredData[startSample];

  rmsPower := 0.0;
  for i := startSample to endSample do
    rmsPower += sqr(globalData^.filteredData[i]);
  rmsPower := sqrt(rmsPower / (endSample - startSample + 1)) - 1.0 / Low(SmallInt);
  rmsPower := 20.0 * log10(rmsPower);

  ComputeBitRange;

  chunks := TChunkList.Create;
  reducedChunks := TChunkList.Create;
end;

destructor TBand.Destroy;
begin
  chunks.Free;
  reducedChunks.Free;

  inherited Destroy;
end;

procedure TBand.ComputeBitRange;
var
  i, hiSmp: Integer;
begin
  hiSmp := 0;
  for i := 0 to frame.encoder.frameSampleCount - 1 do
    hiSmp := max(hiSmp, abs(TEncoder.make16BitSample(srcData[i])));
  bitRange := EnsureRange(1 + ceil(log2(hiSmp + 1.0)), frame.encoder.OutputBitDepth, 16);

  dstBitShift := bitRange - frame.encoder.OutputBitDepth;
end;

procedure TBand.MakeChunks;
var
  i: Integer;
  chunk: TChunk;
begin
  chunks.Capacity := globalData^.chunkCount;
  for i := 0 to globalData^.chunkCount - 1 do
  begin
    chunk := TChunk.Create(Self, i);
    chunk.MakeDstData;
    if not frame.encoder.CRSearch then
      chunk.ComputeFFT(True);
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
        reducedChunk.AddToMix(chunks[i].fft);
        chunks[i].reducedChunk := reducedChunk;
      end;

    reducedChunk.FinalizeMix(True);
  end;

var
  v1, v2, continuityFactor, sl, sh, fdl, fdh: Double;
  i, j, prec: Integer;
  chunk: TChunk;
  XY: TReal2DArray;
  C: TReal2DArray;
begin
  prec := frame.encoder.Precision;
  if (prec = 0) or (desiredChunkCount = globalData^.chunkCount) then Exit;

  //WriteLn('KMeansReduce #', index, ' ', globalData^.desiredChunkCount);

  SetLength(XY, chunks.Count, globalData^.chunkSize * 2 + 4);
  continuityFactor := globalData^.chunkSize / globalData^.chunkSizeUnMin;

  for i := 0 to chunks.Count - 1 do
  begin
    // add first and last sample to features to allow continuity between chunks
    sl := chunks[i].srcData[0];
    sh := chunks[i].srcData[globalData^.chunkSize - 1];
    // add approximate derivatives of first and last sample to features to allow continuity between chunks
    fdl :=  -chunks[i].srcData[0] + chunks[i].srcData[1];
    fdh :=  chunks[i].srcData[globalData^.chunkSize - 1] - chunks[i].srcData[globalData^.chunkSize - 2];

    XY[i, 0] := sl * continuityFactor;
    XY[i, 1] := sh * continuityFactor;
    XY[i, 2] := fdl * continuityFactor;
    XY[i, 3] := fdh * continuityFactor;

    for j := 0 to globalData^.chunkSize - 1 do
    begin
      v1 := chunks[i].fft[j].X;
      v2 := chunks[i].fft[j].Y;

      XY[i, j * 2 + 4] := v1;
      XY[i, j * 2 + 5] := v2;
    end;
  end;

  if frame.encoder.UsePython then
    DoExternalKMeans(XY, desiredChunkCount, 1, prec, False, XYC)
  else
    KMeansGenerate(XY, Length(XY), Length(XY[0]), desiredChunkCount, prec, i, C, XYC);

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

  SetLength(dstData, frame.sampleCount);
  FillQWord(dstData[0], frame.sampleCount, 0);

  pos := 0;
  for i := 0 to chunks.Count - 1 do
  begin
    chunk := chunks[i];

    Assert((frame.encoder.Precision = 0) or (globalData^.chunkCount = desiredChunkCount) or (chunk.reducedChunk <> chunk));

    for j := 0 to globalData^.chunkSize - 1 do
    begin
      smp := TEncoder.makeFloatSample(chunk.reducedChunk.dstData[j], frame.encoder.OutputBitDepth, dstBitShift);
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
var
  j, k, pos, n: Integer;
  acc: Double;
begin
  index := idx;
  band := bnd;

  SetLength(srcData, band.globalData^.chunkSize);

  for j := 0 to band.globalData^.chunkSize - 1 do
  begin
    pos := (idx * band.globalData^.chunkSize + j) * band.globalData^.underSample;

    acc := 0.0;
    n := 0;
    for k := 0 to band.globalData^.underSample - 1 do
    begin
      if pos + k >= band.frame.sampleCount then
        Break;
      acc += band.srcData[pos + k];
      Inc(n);
    end;

    if n = 0 then
      srcData[j] := 0
    else
      srcData[j] := acc / n;
  end;

  reducedChunk := Self;

  mixList := TComplex1DArrayList.Create;
end;

destructor TChunk.Destroy;
begin
  mixList.Free;

  inherited Destroy;
end;

procedure TChunk.ComputeFFT(FromDstData: Boolean);
var
  i: Integer;
  dstf: TReal1DArray;
begin
  if FromDstData then
  begin
    SetLength(dstf, band.globalData^.chunkSize);
    for i := 0 to High(dstf) do
      dstf[i] := (Integer(dstData[i]) + Low(ShortInt)) / High(ShortInt);
    FFTR1D(dstf, band.globalData^.chunkSize, fft);
  end
  else
  begin
    FFTR1D(srcData, band.globalData^.chunkSize, fft);
  end;
end;

procedure TChunk.InitMix;
begin
  mixList.Clear;
end;

procedure TChunk.AddToMix(afft: TComplex1DArray);
begin
  mixList.Add(afft);
end;

function CompareComplexX(const d1,d2): integer;
var
  c1 : Complex absolute d1;
  c2 : Complex absolute d2;
begin
  Result := CompareValue(c1.X, c2.X);
end;

procedure TChunk.FinalizeMix(IntoDstData: Boolean);
var
  i, k: Integer;
  acca: TComplex1DArray;
  acc: Complex;
  dstf: TReal1DArray;
begin
  if not Assigned(mixList) or (mixList.Count = 0) then
    Exit;

  SetLength(fft, band.globalData^.chunkSize);
  SetLength(acca, mixList.Count);

  for k := 0 to band.globalData^.chunkSize - 1 do
  begin

    acc.X := 0;
    acc.Y := 0;
    if mixList.Count >= 7 then
    begin
      // remove outliers

      for i := 0 to mixList.Count - 1 do
        acca[i]:= mixList[i][k];

      anysort.AnySort(acca[0], mixList.Count, SizeOf(Complex), @CompareComplexX);

      for i := 1 to mixList.Count - 2 do
      begin
        acc.X += acca[i].X;
        acc.Y += acca[i].Y;
      end;

      acc.X /= mixList.Count - 2;
      acc.Y /= mixList.Count - 2;
    end
    else
    begin
      for i := 0 to mixList.Count - 1 do
      begin
        acc.X += mixList[i][k].X;
        acc.Y += mixList[i][k].Y;
      end;

      acc.X /= mixList.Count;
      acc.Y /= mixList.Count;
    end;

    fft[k] := acc;
  end;

  if IntoDstData then
  begin
    FFTR1DInv(fft, band.globalData^.chunkSize, dstf);
    SetLength(dstData, band.globalData^.chunkSize);
    for i := 0 to band.globalData^.chunkSize - 1 do
      dstData[i] := EnsureRange(round(dstf[i] * High(ShortInt)), Low(ShortInt), High(ShortInt)) - Low(ShortInt);
  end
  else
  begin
    FFTR1DInv(fft, band.globalData^.chunkSize, srcData);
    MakeDstData;
  end;
end;

procedure TChunk.MakeDstData;
var
  i: Integer;
begin
  SetLength(dstData, band.globalData^.chunkSize);
  for i := 0 to band.globalData^.chunkSize - 1 do
    dstData[i] := TEncoder.makeOutputSample(srcData[i], band.frame.encoder.OutputBitDepth, band.dstBitShift);
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
    sampleCount := (fs.Size - fs.Position) div SizeOf(SmallInt);
    SetLength(srcData, sampleCount);
    FillWord(srcData[0], sampleCount, 0);
    for i := 0 to sampleCount - 1 do
      srcData[i] := SmallInt(fs.ReadWord);
  finally
    fs.Free;
  end;

  sampleRate := PInteger(@srcHeader[$18])^;
  writeln('sampleRate = ', sampleRate);
  WriteLn('srcDataCount = ', sampleCount);
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
    fs.WriteBuffer(dstData[0], sampleCount * 2);
  finally
    fs.Free;
  end;
end;

procedure TEncoder.SaveRSC;
var
  i: Integer;
  fs: TFileStream;
  best, cur: TMemoryStream;
  fn: String;
begin
  fs := nil;
  fn := ChangeFileExt(outputFN, '.rsc');
  best := TMemoryStream.Create;
  cur := TMemoryStream.Create;
  try
    WriteLn('Save ', fn);

    for i := 0 to 1000*0 do
    begin
      cur.Clear;
      SaveStream(cur);
      cur.Position := 0;

      if (best.Size = 0) or (cur.Size < best.Size) then
      begin
        cur.Position := 0;
        best.Clear;
        best.CopyFrom(cur, cur.Size);
        WriteLn('Iteration #', i, ', best: ', best.Size);
      end;
    end;

    fs := TFileStream.Create(fn, fmCreate);
    best.Position := 0;
    fs.CopyFrom(cur, fs.Size);
  finally
    fs.Free;
    best.Free;
    cur.Free;
  end;
end;

procedure TEncoder.SaveStream(AStream: TStream);
const
  CacheSize = 256;
var
  blockSize, blockCount: Integer;
  chunksPerBlock: array[0 .. BandCount - 1] of Integer;
  cache: array[0 .. CacheSize - 1, 0..2] of Integer;
  blockCacheIdxs: array[0 .. BandCount - 1, 0 .. CacheSize - 1] of Integer;

  function IsCacheIdxUsed(idx: Integer): Boolean;
  var
    i, j: Integer;
  begin
    Result := idx = 0; // cache idx #0 doesn't exist (play command)
    for i := 0 to BandCount - 1 do
      for j := 0 to chunksPerBlock[i] - 1 do
        if blockCacheIdxs[i, j] = idx then
          Exit(True);
  end;

  function AllocCacheIdx(chunk: TChunk): Integer;
  var
    i, ratio, lru: Integer;
  begin
    ratio := chunk.band.globalData^.chunkSize div MinChunkSize;

{$if false}
    repeat
      Result := random(CacheSize div ratio) * ratio;
    until not IsCacheIdxUsed(Result); // don't use a cache index that is needed for this block
{$else}
    Result := CacheSize - 1;
    lru := MaxInt;
    for i := 0 to CacheSize div ratio - 1 do
      if cache[i * ratio, 2] < lru then
      begin
        Result := i * ratio;
        lru := cache[Result, 2];
      end;
{$endif}

    Assert(not IsCacheIdxUsed(Result));

    cache[Result, 0] := chunk.band.index;
    cache[Result, 1] := chunk.index;
    for i := 1 to ratio - 1 do
    begin
      cache[Result + i, 0] := -chunk.band.index;
      cache[Result + i, 1] := -chunk.index;
    end;

    //AStream.WriteWord($aa55);
    AStream.WriteByte(Result);
    AStream.Write(chunk.dstData[0], chunk.band.globalData^.chunkSize);
  end;

  function GetOrAllocCacheIdx(chunk: TChunk; lru: Integer): Integer;
  var
    i: Integer;
  begin
    Result := -1;
    for i := 0 to CacheSize - 1 do
      if (cache[i, 0] = chunk.band.index) and (cache[i, 1] = chunk.index) then
      begin
        Result := i;
        Break;
      end;

    if Result < 0 then
      Result := AllocCacheIdx(chunk);

    cache[Result, 2] := lru;
  end;

//var
//  i, j, k: Integer;
  //pp: Integer;
begin
  //blockSize := bands[0].chunkSize * bands[0].underSample;
  //blockCount := (sampleCount - 1) div blockSize + 1;
  //
  //for i := 0 to BandCount - 1 do
  //  chunksPerBlock[i] := blockSize div (bands[i].chunkSize * bands[i].underSample);
  //
  //FillChar(cache[0, 0], SizeOf(cache), Byte(-1));
  //cache[0, 2] := MaxInt; // cache idx #0 doesn't exist (play command)
  //
  ////pp := AStream.Position;
  //for i := 0 to blockCount - 1 do
  //begin
  //  FillChar(blockCacheIdxs[0, 0], SizeOf(blockCacheIdxs), Byte(-1));
  //
  //  for j := 0 to BandCount - 1 do
  //    for k := 0 to chunksPerBlock[j] - 1 do
  //      blockCacheIdxs[j, k] := GetOrAllocCacheIdx(bands[j].chunkList[Min(bands[j].chunkCount - 1, i * chunksPerBlock[j] + k)].reducedChunk, i);
  //
  //  AStream.WriteByte(0); // play command
  //  for j := 0 to BandCount - 1 do
  //    for k := 0 to chunksPerBlock[j] - 1 do
  //      AStream.WriteByte(blockCacheIdxs[j, k]);
  //
  //  //writeln(AStream.Position - pp);
  //  //pp := AStream.Position;
  //end;
end;

procedure TEncoder.SaveBandWAV(index: Integer; fn: String);
var
  fs: TFileStream;
begin
  WriteLn('SaveBandWAV #', index, ' ', fn);

  fs := TFileStream.Create(fn, fmCreate or fmShareDenyWrite);
  try
    fs.WriteBuffer(srcHeader[0], SizeOf(srcHeader));
    fs.WriteBuffer(bandData[index].dstData[0], sampleCount * 2);
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

    hc := min(HighCut, sampleRate / 2);
    ratioP := log2(LowCut / hc) / BandCount + 0.5;
    ratioL := hc / sampleRate;

    if i = 0 then
      bnd.fcl := LowCut / sampleRate
    else
      bnd.fcl := power(2.0, (BandCount - i) * ratioP) * ratioL;

    if i = BandCount - 1 then
      bnd.fch := hc / sampleRate
    else
      bnd.fch := power(2.0, (BandCount - 1 - i) * ratioP) * ratioL;

    // undersample if the band high freq is a lot lower than nyquist

    bnd.chunkSize := round(intpower(2.0, floor(-log2(bnd.fcl))));
    bnd.underSample := round(intpower(2.0, floor(-log2(bnd.fch)) - 2));
    bnd.underSample := Max(1, bnd.underSample);
    bnd.chunkSize := bnd.chunkSize div bnd.underSample;

    bnd.underSampleUnMin := bnd.underSample;
    bnd.chunkSizeUnMin := bnd.chunkSize;
    if bnd.chunkSize < minChunkSize then
    begin
      bnd.underSample := max(1, (bnd.underSample * bnd.chunkSize) div minChunkSize);
      bnd.chunkSize := minChunkSize;
    end;

    bandData[i] := bnd;
  end;
end;

procedure TEncoder.MakeBandSrcData(AIndex: Integer);
var
  j: Integer;
  bnd: TBandGlobalData;
  cic: TSecondOrderCIC;
  cicSmp, cr1l, cr2l: Double;
begin
  bnd := bandData[AIndex];

  SetLength(bnd.filteredData, sampleCount);
  for j := 0 to sampleCount - 1 do
    bnd.filteredData[j] := makeFloatSample(srcData[j]);

  // band pass the samples
  bnd.filteredData := DoBPFilter(bnd.fcl, bnd.fch, BandTransFactor, bnd.filteredData);
  SetLength(bnd.filteredData, frameSampleCount * frameCount);

  for j := sampleCount to frameSampleCount * frameCount - 1 do
    bnd.filteredData[j] := 0;

  if bnd.underSample > 1 then
  begin
    // compensate for decoder altering the pass band
    cic := TSecondOrderCIC.Create(bnd.underSample * 2, BandDealiasSecondOrder);

    cr1l := CICCorrRatio[BandDealiasSecondOrder, 0];
    cr2l := CICCorrRatio[BandDealiasSecondOrder, 1];
    if CRSearch then
    begin
      cr1l := CR1;
      cr2l := CR2;
    end;

    try
      for j := -cic.Ratio + 1 to High(bnd.filteredData) do
      begin
        cicSmp := cic.ProcessSample(bnd.filteredData[Min(High(bnd.filteredData), j + cic.Ratio - 1)]);
        if j >= 0 then
          bnd.filteredData[j] += bnd.filteredData[j] * cr1l - cicSmp * cr2l;
      end;
    finally
      cic.Free;
    end;
  end;

  bandData[AIndex] := bnd;
end;

procedure TEncoder.MakeBandGlobalData2;
var
  wgt, fsq: Double;
  i: Integer;
  bnd: TBandGlobalData;
begin
  for i := 0 to BandCount - 1 do
  begin
    bnd := bandData[i];

    bnd.chunkCount := (frameSampleCount - 1) div (bnd.chunkSize * bnd.underSample) + 1;

    fsq := bnd.fcl * bnd.fch * sampleRate * sampleRate;

    if BandBWeighting then
      // B-weighting
      wgt := sqr(12194.0) * sqrt(fsq) * fsq / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)) * sqrt(fsq + sqr(158.5)))
    else
      // A-weighting
      wgt := sqr(12194.0) * sqr(fsq) / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)) * sqrt((fsq + sqr(107.7)) * (fsq + sqr(737.9))));

    bnd.weight := power(wgt, BandWeightingApplyPower);

    bandData[i] := bnd;
  end;
end;


procedure TEncoder.MakeFrames;

  procedure DoBand(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  begin
    MakeBandSrcData(AIndex);
  end;

  procedure DoFrame(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    frm: TFrame;
  begin
    frm := frames[AIndex];

    frm.MakeBands;
  end;

var
  i, blockSampleCount, blockChunkCount, blockCount, curStart, curEnd: Integer;
  frm: TFrame;
begin
  MakeBandGlobalData;

  blockSampleCount := 0;
  for i := 0 to BandCount - 1 do
    if bandData[i].underSample > blockSampleCount then
      blockSampleCount := bandData[i].underSample * bandData[i].chunkSize;

  blockChunkCount := 0;
  for i := 0 to BandCount - 1 do
    blockChunkCount += blockSampleCount div (bandData[i].underSample * bandData[i].chunkSize);

  blockCount := (sampleCount - 1) div blockSampleCount + 1;

  // pass 1
  frameByteSize := ChunksPerFrame * MinChunkSize;
  projectedByteSize := ceil((sampleCount / sampleRate) * BitRate * (1024 / 8));
  frameCount :=  (projectedByteSize - blockCount * blockChunkCount - 1) div frameByteSize + 1;
  frameSampleCount := (sampleCount - 1) div frameCount + 1;
  frameSampleCount := ((frameSampleCount - 1) div blockSampleCount + 1) * blockSampleCount;

  // pass 2
  frameByteSize += frameSampleCount div blockSampleCount * blockChunkCount;
  projectedByteSize := frameByteSize * frameCount;

  MakeBandGlobalData2;

  if verbose then
  begin
    writeln('frameByteSize = ', frameByteSize);
    writeln('frameCount = ', frameCount);
    writeln('frameSampleCount = ', frameSampleCount);
    writeln('projectedByteSize = ', projectedByteSize);
  end;

  for i := 0 to BandCount - 1 do
     WriteLn('Band #', i, ' (', round(bandData[i].fcl * sampleRate), ' Hz .. ', round(bandData[i].fch * sampleRate), ' Hz); ', bandData[i].chunkSize, ' * ', bandData[i].chunkCount, '; ', bandData[i].underSample);

  ProcThreadPool.DoParallelLocalProc(@DoBand, 0, BandCount - 1, nil);

  for i := 0 to frameCount - 1 do
  begin
    curStart := i * frameSampleCount;
    curEnd := min(sampleCount, (i + 1) * frameSampleCount) - 1;

    frm := TFrame.Create(Self, i, curStart, curEnd);
    frames.Add(frm);
  end;

  ProcThreadPool.DoParallelLocalProc(@DoFrame, 0, frameCount - 1, nil);
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
    coeffs := DoFilterCoeffs(fcl, transFactor * (exp(fcl) - 1.0), True);
    Result := DoFilter(Result, coeffs);
  end;

  if fch < 0.5 then
  begin
    coeffs := DoFilterCoeffs(fch, transFactor * (exp(fch) - 1.0), False);
    Result := DoFilter(Result, coeffs);
  end;
end;

constructor TEncoder.Create(InFN, OutFN: String);
begin
  inputFN := InFN;
  outputFN := OutFN;

  BitRate := 128;
  Precision := 10;
  BandTransFactor := 0.1;
  LowCut := 32.0;
  HighCut := 16000.0;
  BandDealiasSecondOrder := True;
  BandBWeighting := True;
  BandWeightingApplyPower := 0.0;
  OutputBitDepth := 8;
  MinChunkSize := 16;
  ChunksPerFrame:= 256;
  UsePython := False;

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
  cic: array[0 .. BandCount - 1] of TSecondOrderCIC;
  offset: array[0 .. BandCount - 1] of Integer;
begin
  WriteLn('MakeDstData');

  SetLength(dstData, frameSampleCount * frames.Count);
  FillWord(dstData[0], Length(dstData), 0);

  for i := 0 to BandCount - 1 do
  begin
    SetLength(bandData[i].dstData, Length(dstData));
    FillWord(bandData[i].dstData[0], Length(dstData), 0);

    if bandData[i].underSample > 1 then
    begin
      cic[i] := TSecondOrderCIC.Create(bandData[i].underSample * 2, BandDealiasSecondOrder);
      offset[i] := -cic[i].Ratio + 1;
    end
    else
    begin
      cic[i] := nil;
      offset[i] := 0;
    end;
  end;

  pos := 0;
  try
    for k := 0 to frames.Count - 1 do
      for i := 0 to frames[k].sampleCount - 1 do
      begin
        for j := 0 to BandCount - 1 do
        begin
          smp := frames[k].bands[j].dstData[i];

          if Assigned(cic[j]) then
            smp := cic[j].ProcessSample(smp);

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
      cic[i].Free;
  end;
end;

procedure TEncoder.SearchBestCorrRatios;

  function DoOne(ACR1, ACR2: Double; Verbose: Boolean): Double;
  begin
    CR1 := CICCorrRatio[BandDealiasSecondOrder, 0];
    CR2 := ACR1;

    MakeFrames;
    MakeDstData;
    SaveWAV;
    Result := DoExternalEAQUAL(ParamStr(1), ParamStr(2), False, False, 2048);

    if Verbose then
      Write(#13'(', FormatFloat('0.00000', CR1), ', ', FormatFloat('0.00000', CR2), ')'#9, FormatFloat('0.00000', Result), '                ');
  end;

var
  N : AlglibInteger;
  State : MinLBFGSState;
  Rep : MinLBFGSReport;
  S : TReal1DArray;
  H : Double;
begin
  CRSearch := True;

  H := ifthen(BandDealiasSecondOrder, 1e-4, 1e-4);
  N := 1;
  SetLength(S, N);
  S[0] := ifthen(BandDealiasSecondOrder, 1.8, 0.5);
  MinLBFGSCreate(N, N, S, State);
  MinLBFGSSetCond(State, H, 0.0, 0.0, 0);
  while MinLBFGSIteration(State) do
  begin
    if State.NeedFG then
    begin
      State.F := DoOne(State.X[0], State.X[0], True);

{$if false}
      State.G[0] := (
          -DoOne(State.X[0] + 2 * H, State.X[1], False) +
          8 * DoOne(State.X[0] + H, State.X[1], False) +
          -8 * DoOne(State.X[0] - H, State.X[1], False) +
          DoOne(State.X[0] - 2 * H, State.X[1], False)
        ) / (12 * H);
{$else}
      State.G[0] := (
          DoOne(State.X[0] + H, State.X[1], False) -
          DoOne(State.X[0] - H, State.X[1], False)
        ) / (2 * H);
{$endif}
    end;
  end;
  MinLBFGSResults(State, S, Rep);

  WriteLn;
  WriteLn('Best: ');
  DoOne(S[0], S[1], True);
end;

function TEncoder.DoFilterCoeffs(fc, transFactor: Double; HighPass: Boolean): TDoubleDynArray;
var
  sinc, win, sum: Double;
  i, N: Integer;
begin
  N := ceil(4.6 / transFactor);
  if (N mod 2) = 0 then N += 1;

  //writeln('DoFilterCoeffs ', ifthen(HighPass, 'HP', 'LP'), ' ', FloatToStr(sampleRate * fc), ' ', N);

  SetLength(Result, N);
  sum := 0;
  for i := 0 to N - 1 do
  begin
    sinc := 2.0 * fc * (i - (N - 1) / 2.0) * pi;
    if IsZero(sinc) then
      sinc := 1.0
    else
      sinc := sin(sinc) / sinc;

    // blackman window
    win := 0.42 - 0.5 * cos(2 * pi * i / (N - 1)) + 0.08 * cos(4 * pi * i / (N - 1));

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
  smp16 := make16BitSample(smp);
  smp16 := (smp16 - (sign(smp16) shl (max(0, bitShift - 1)))) div (1 shl bitShift);
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
  smp16 := smp16 * (1 shl bitShift);
  Result := makeFloatSample(smp16);
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

function TEncoder.ComputeEAQUAL(chunkSz: Integer; UseDIX: Boolean; const smpRef, smpTst: TDoubleDynArray): Double;
var
  i: Integer;
  FNRef, FNTst: String;
  ms: TMemoryStream;
begin

  FNRef := GetTempFileName('', 'ref-'+IntToStr(GetCurrentThreadId)+'.wav');
  FNTst := GetTempFileName('', 'tst-'+IntToStr(GetCurrentThreadId)+'.wav');

  ms := TMemoryStream.Create;
  try
    ms.Write(srcHeader[0], $28);
    ms.WriteDWord(chunkSz * SizeOf(SmallInt));

    for i := 0 to chunkSz - 1 do
      ms.WriteWord(make16BitSample(smpRef[i]) - Low(SmallInt));

    ms.SaveToFile(FNRef);
    ms.Clear;

    ms.Write(srcHeader[0], $28);
    ms.WriteDWord(chunkSz * SizeOf(SmallInt));

    for i := 0 to chunkSz - 1 do
      ms.WriteWord(make16BitSample(smpTst[i]) - Low(SmallInt));

    ms.SaveToFile(FNTst);
  finally
    ms.Free;
  end;

  Result := DoExternalEAQUAL(FNRef, FNTst, False, UseDIX, chunkSz * 2);

  DeleteFile(FNRef);
  DeleteFile(FNTst);
end;

function TEncoder.ComputeEAQUALMulti(chunkSz: Integer; UseDIX: Boolean; const smpRef: TDoubleDynArray;
  smpTst: TDoubleDynArray2): TDoubleDynArray;
var
  i, j: Integer;
  FNRef, FNTst: String;
  zeroes: TSmallIntDynArray;
  ms: TMemoryStream;
begin
  if Length(smpTst) = 0 then
    Exit(nil);

  SetLength(zeroes, chunkSz);
  FillWord(zeroes[0], chunkSz, 0);

  FNRef := GetTempFileName('', 'ref-'+IntToStr(GetCurrentThreadId)+'.wav');
  FNTst := GetTempFileName('', 'tst-'+IntToStr(GetCurrentThreadId)+'.wav');

  ms := TMemoryStream.Create;
  try
    ms.Write(srcHeader[0], $28);
    ms.WriteDWord(chunkSz * 2 * SizeOf(SmallInt) * Length(smpTst));
    for i := 0 to High(smpTst) do
    begin
      for j := 0 to chunkSz - 1 do
        ms.WriteWord(make16BitSample(smpRef[j]) - Low(SmallInt));
      ms.Write(zeroes[0], chunkSz * SizeOf(SmallInt));
    end;

    ms.SaveToFile(FNRef);
    ms.Clear;

    ms.Write(srcHeader[0], $28);
    ms.WriteDWord(chunkSz * 2 * SizeOf(SmallInt) * Length(smpTst));
    for i := 0 to High(smpTst) do
    begin
      for j := 0 to chunkSz - 1 do
        ms.WriteWord(make16BitSample(smpTst[i, j]) - Low(SmallInt));
      ms.Write(zeroes[0], chunkSz * SizeOf(SmallInt));
    end;

    ms.SaveToFile(FNTst);
  finally
    ms.Free;
  end;

  Result := DoExternalEAQUALMulti(FNRef, FNTst, UseDIX, Length(smpTst), chunkSz * 2);

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
      WriteLn(#9'-bwa'#9'use A-weighting for bands (otherwise B-weighting)');
      WriteLn(#9'-bwp'#9'band weighting apply power');
      WriteLn(#9'-obd'#9'output bit depth (1-8)');
      WriteLn(#9'-mcs'#9'minimum chunk size');
      WriteLn(#9'-v'#9'verbose K-means');
      WriteLn(#9'-py'#9'use Python script for clustering');
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
      enc.BandBWeighting :=  ParamStart('-bwa') = -1;
      enc.BandWeightingApplyPower := round(ParamValue('-bwp', enc.BandWeightingApplyPower));
      enc.OutputBitDepth :=  round(ParamValue('-obd', enc.OutputBitDepth));
      enc.MinChunkSize :=  round(ParamValue('-mcs', enc.MinChunkSize));
      enc.verbose := ParamStart('-v') <> -1;
      enc.UsePython := ParamStart('-py') <> -1;

      WriteLn('BitRate = ', FloatToStr(enc.BitRate));
      WriteLn('Precision = ', enc.Precision);
      WriteLn('LowCut = ', FloatToStr(enc.LowCut));
      WriteLn('HighCut = ', FloatToStr(enc.HighCut));
      WriteLn('BandTransFactor = ', FloatToStr(enc.BandTransFactor));
      WriteLn('BandDealiasSecondOrder = ', BoolToStr(enc.BandDealiasSecondOrder, True));
      WriteLn('BandBWeighting = ', BoolToStr(enc.BandBWeighting, True));
      WriteLn('BandWeightingApplyPower = ', FloatToStr(enc.BandWeightingApplyPower));
      WriteLn('OutputBitDepth = ', enc.OutputBitDepth);
      WriteLn('MinChunkSize = ', enc.MinChunkSize);
      WriteLn('UsePython = ', BoolToStr(enc.UsePython, True));
      WriteLn;

      enc.Load;
      enc.MakeFrames;
      enc.MakeDstData;
      enc.SaveWAV;
      for i := 0 to BandCount - 1 do
        enc.SaveBandWAV(i, ChangeFileExt(enc.outputFN, '-' + IntToStr(i) + '.wav'));
      enc.SaveRSC;

      DoExternalEAQUAL(ParamStr(1), ParamStr(2), True, True, 2048);

      if ParamStart('-scr') <> -1 then
        enc.SearchBestCorrRatios;

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

