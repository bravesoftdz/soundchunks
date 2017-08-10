program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, yakmo, ap, fft, conv, anysort, minlbfgs, kmeans, correlation;

const
  BandCount = 4;
  InputSNRDb = -90.3; // 16bit
  CICCorrRatio: array[Boolean{2nd order?}, 0..1] of Double = (
    (0.503, 0.816),
    (1.006, 0.949)
  );

type
  TEncoder = class;
  TFrame = class;
  TBand = class;
  TChunk = class;

  TBandGlobalData = record
    fcl, fch: Double;
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
    R, R2, pos, pos2: Integer;
    smps, vs: TDoubleDynArray;
  public
    Ratio: Integer;
    SecondOrder: Boolean;

    constructor Create(ARatio: Integer; ASecondOrder: Boolean);

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
    tmpIndex: Integer;

    mixList: TChunkList;

    origSrcData: PDouble;
    srcData: TDoubleDynArray;
    fft: TComplex1DArray;
    dstData: TByteDynArray;

    constructor Create(bnd: TBand; idx: Integer);
    destructor Destroy; override;

    procedure ComputeFFT(FromDstData: Boolean);
    procedure InitMix;
    procedure AddToMix(chk: TChunk);
    procedure FinalizeMix(IntoDstData: Boolean);
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
    OutputBitDepth: Integer; // max 8Bits
    MinChunkSize: Integer;
    MaxFrameSize: Integer;
    MaxChunksPerBand: Integer;
    UsePython: Boolean;

    SampleRate: Integer;
    SampleCount: Integer;
    projectedByteSize, frameCount, frameSampleCount: Integer;
    verbose: Boolean;

    BWSearch, CRSearch: Boolean;
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

    procedure PrepareFrames;
    procedure MakeFrames;
    procedure MakeDstData;

    procedure SearchBestBandWeighting;
    procedure SearchBestCorrRatios;

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

constructor TFrame.Create(enc: TEncoder; idx: Integer; startSample, endSample: Integer);
var
  i: Integer;
begin
  encoder := enc;
  index := idx;

  BandRMSWeighting := True;
  BandBWeighting := True;
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
  i, allSz: Integer;
  sz, dbh: Double;
  wgt, fsq: Double;
  full, compressed: Boolean;
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


  allSz := 0;
  sz := 1;
  repeat
    FrameSize := allSz;
    for i := 0 to BandCount - 1 do
      bands[i].desiredChunkCount := dcc[i];

    full := True;
    allSz := 0;
    for i := 0 to BandCount - 1 do
    begin
      bnd := bands[i];

      dcc[i] := round(sz * bnd.weight);
      dcc[i] := EnsureRange(dcc[i], 1, min(encoder.MaxChunksPerBand, bnd.globalData^.chunkCount));

      compressed := dcc[i] < min(encoder.MaxChunksPerBand, bnd.globalData^.chunkCount);
      full := full and not compressed;

      allSz += dcc[i] * bnd.globalData^.chunkSize + Ord(compressed) * bnd.globalData^.chunkCount + SizeOf(Byte);
    end;
    sz += 0.1;

  until (allSz > encoder.MaxFrameSize) or full;


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

  if not encoder.verbose or encoder.BWSearch then
    Write('.');
end;

constructor TSecondOrderCIC.Create(ARatio: Integer; ASecondOrder: Boolean);
begin
  SecondOrder := ASecondOrder;
  Ratio := ARatio;
  R := max(1, ARatio);
  R2 := R div 2;
  SetLength(smps, R);
  SetLength(vs, R2);
  v := 0;
  v2 := 0;
  pos := 0;
  pos2 := 0;
end;

function TSecondOrderCIC.ProcessSample(Smp: Double): Double;
begin
  Smp /= R;

  v := v + smp - smps[pos];
  smps[pos] := smp;
  pos := (pos + 1) mod R;

  vv := v / R2;
  v2 := v2 + vv - vs[pos2];
  vs[pos2] := vv;
  pos2 := (pos2 + 1) mod R2;


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
      chunk.ComputeFFT(False);
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
    fdl :=  -1.5 * chunks[i].srcData[0] +
            2.0 * chunks[i].srcData[1] +
            -0.5 * chunks[i].srcData[2];
    fdh :=  1.5 * chunks[i].srcData[globalData^.chunkSize - 1] +
            -2.0 * chunks[i].srcData[globalData^.chunkSize - 2] +
            0.5 * chunks[i].srcData[globalData^.chunkSize - 3];

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
      smp := TEncoder.makeFloatSample(chunk.reducedChunk.dstData[j], frame.encoder.OutputBitDepth, chunk.dstBitShift);
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

procedure TChunk.ComputeFFT(FromDstData: Boolean);
var
  i: Integer;
  ffb: TReal1DArray;
begin
  SetLength(ffb, band.globalData^.chunkSize);
  if FromDstData then
  begin
    for i := 0 to High(ffb) do
      ffb[i] := (Integer(dstData[i]) + Low(ShortInt)) / High(ShortInt);
  end
  else
  begin
    for i := 0 to High(ffb) do
      ffb[i] := srcData[i] * (1 shl (16 - bitRange));
  end;
  FFTR1D(ffb, band.globalData^.chunkSize, fft);
end;

procedure TChunk.InitMix;
begin
  mixList.Clear;
end;

procedure TChunk.AddToMix(chk: TChunk);
begin
  mixList.Add(chk);
end;

procedure TChunk.FinalizeMix(IntoDstData: Boolean);
var
  i, k, sz: Integer;
  acc: Integer;
  accf: Double;
  osrc: TDoubleDynArray;
begin
  if not Assigned(mixList) or (mixList.Count = 0) then
    Exit;

  if IntoDstData then
  begin
    SetLength(dstData, band.globalData^.chunkSize);

    for k := 0 to band.globalData^.chunkSize - 1 do
    begin
      acc := 0;

      for i := 0 to mixList.Count - 1 do
        acc += mixList[i].dstData[k];

      dstData[k] := acc div mixList.Count;
    end;
  end
  else
  begin
    sz := band.globalData^.chunkSize * band.globalData^.underSample * 2;
    SetLength(osrc, sz * 3);

    for k := 0 to High(osrc) do
    begin
      accf := 0.0;

      for i := 0 to mixList.Count - 1 do
        accf += mixList[i].origSrcData[k] * (1 shl (16 - mixList[i].bitRange));

      osrc[k] := accf / mixList.Count;
    end;

    MakeSrcData(@osrc[0]);

    bitRange := 16;
    dstBitShift := 8;

    MakeDstData;
  end;
end;

procedure TChunk.ComputeBitRange;
var
  i, hiSmp: Integer;
begin
  hiSmp := 0;
  for i := 0 to band.globalData^.chunkSize - 1 do
    hiSmp := max(hiSmp, abs(TEncoder.make16BitSample(srcData[i])));
  bitRange := EnsureRange(1 + ceil(log2(hiSmp + 1.0)), band.frame.encoder.OutputBitDepth, 16);
  dstBitShift := bitRange - band.frame.encoder.OutputBitDepth;
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
    dstData[i] := TEncoder.makeOutputSample(srcData[i], band.frame.encoder.OutputBitDepth, dstBitShift);
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
//const
//  CacheSize = 256;
//var
//  blockSize, blockCount: Integer;
//  chunksPerBlock: array[0 .. BandCount - 1] of Integer;
//  cache: array[0 .. CacheSize - 1, 0..2] of Integer;
//  blockCacheIdxs: array[0 .. BandCount - 1, 0 .. CacheSize - 1] of Integer;
//
//  function IsCacheIdxUsed(idx: Integer): Boolean;
//  var
//    i, j: Integer;
//  begin
//    Result := idx = 0; // cache idx #0 doesn't exist (play command)
//    for i := 0 to BandCount - 1 do
//      for j := 0 to chunksPerBlock[i] - 1 do
//        if blockCacheIdxs[i, j] = idx then
//          Exit(True);
//  end;
//
//  function AllocCacheIdx(chunk: TChunk): Integer;
//  var
//    i, ratio, lru: Integer;
//  begin
//    ratio := chunk.band.globalData^.chunkSize div MinChunkSize;
//
//{$if false}
//    repeat
//      Result := random(CacheSize div ratio) * ratio;
//    until not IsCacheIdxUsed(Result); // don't use a cache index that is needed for this block
//{$else}
//    Result := CacheSize - 1;
//    lru := MaxInt;
//    for i := 0 to CacheSize div ratio - 1 do
//      if cache[i * ratio, 2] < lru then
//      begin
//        Result := i * ratio;
//        lru := cache[Result, 2];
//      end;
//{$endif}
//
//    Assert(not IsCacheIdxUsed(Result));
//
//    cache[Result, 0] := chunk.band.index;
//    cache[Result, 1] := chunk.index;
//    for i := 1 to ratio - 1 do
//    begin
//      cache[Result + i, 0] := -chunk.band.index;
//      cache[Result + i, 1] := -chunk.index;
//    end;
//
//    //AStream.WriteWord($aa55);
//    AStream.WriteByte(Result);
//    AStream.Write(chunk.dstData[0], chunk.band.globalData^.chunkSize);
//  end;
//
//  function GetOrAllocCacheIdx(chunk: TChunk; lru: Integer): Integer;
//  var
//    i: Integer;
//  begin
//    Result := -1;
//    for i := 0 to CacheSize - 1 do
//      if (cache[i, 0] = chunk.band.index) and (cache[i, 1] = chunk.index) then
//      begin
//        Result := i;
//        Break;
//      end;
//
//    if Result < 0 then
//      Result := AllocCacheIdx(chunk);
//
//    cache[Result, 2] := lru;
//  end;
//
//var
//  i, j, k: Integer;
  //pp: Integer;
begin
  //blockSize := bands[0].chunkSize * bands[0].underSample;
  //blockCount := (SampleCount - 1) div blockSize + 1;
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
    ratioP := log2(LowCut / hc) / BandCount + 0.5;
    ratioL := hc / SampleRate;

    if i = 0 then
      bnd.fcl := LowCut / SampleRate
    else
      bnd.fcl := power(2.0, (BandCount - i) * ratioP) * ratioL;

    if i = BandCount - 1 then
      bnd.fch := hc / SampleRate
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

  SetLength(bnd.filteredData, SampleCount);
  for j := 0 to SampleCount - 1 do
    bnd.filteredData[j] := makeFloatSample(srcData[j]);

  // band pass the samples
  bnd.filteredData := DoBPFilter(bnd.fcl, bnd.fch, BandTransFactor, bnd.filteredData);
  SetLength(bnd.filteredData, frameSampleCount * frameCount);

  for j := SampleCount to frameSampleCount * frameCount - 1 do
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

procedure TEncoder.PrepareFrames;

  procedure DoBand(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  begin
    MakeBandSrcData(AIndex);
  end;

var
  i, blockSampleCount, blockChunkCount, curStart, curEnd: Integer;
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

  // pass 1
  projectedByteSize := ceil((SampleCount / SampleRate) * BitRate * (1024 / 8));
  frameCount :=  (projectedByteSize - 1) div MaxFrameSize + 1;
  frameSampleCount := (SampleCount - 1) div frameCount + 1;
  frameSampleCount := ((frameSampleCount - 1) div blockSampleCount + 1) * blockSampleCount;

  // pass 2
  projectedByteSize := MaxFrameSize * frameCount;

  MakeBandGlobalData2;

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

procedure TEncoder.MakeBandGlobalData2;
var
  i: Integer;
  bnd: TBandGlobalData;
begin
  for i := 0 to BandCount - 1 do
  begin
    bnd := bandData[i];
    bnd.chunkCount := (frameSampleCount - 1) div (bnd.chunkSize * bnd.underSample) + 1;
    bandData[i] := bnd;
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
  BandTransFactor := 0.5;
  LowCut := 32.0;
  HighCut := 18000.0;
  BandDealiasSecondOrder := True;
  OutputBitDepth := 8;
  MinChunkSize := 16;
  MaxFrameSize:= 7168;
  MaxChunksPerBand := 128;
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
  if not (BWSearch or CRSearch) then
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
      for i := 0 to frames[k].SampleCount - 1 do
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

procedure TEncoder.SearchBestBandWeighting;
const
  TestCount = 11;
  bwr: array[0..10] of Boolean = (True, True, True, True, True, True, True, False, False, False, False);
  bwa: array[0..10] of Boolean = (False, False, False, True, True, False, False, False, False, False, True);
  bwp: array[0..10] of Double = (0, 1, -1, 1, -1, 2, -2, 0, 1, 2, 1);

var
  results: array[0..TestCount - 1] of TDoubleDynArray;

  procedure DoEAQ(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    ref, tst: TSmallIntDynArray;
  begin
    ref := Copy(srcData, AIndex * frameSampleCount, frameSampleCount);
    SetLength(ref, frameSampleCount);
    tst := Copy(dstData, AIndex * frameSampleCount, frameSampleCount);
    results[PtrInt(AData), AIndex] := ComputeEAQUAL(frameSampleCount, False, False, ref, tst);
  end;

var
  i, j: Integer;
  best: Double;
begin
  BWSearch := True;
  try
    for i := 0 to TestCount - 1 do
    begin
      PrepareFrames;

      for j := 0 to frameCount - 1 do
      begin
        frames[j].BandRMSWeighting := bwr[i];
        frames[j].BandBWeighting := not bwa[i];
        frames[j].BandWeightingApplyPower := bwp[i];
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
    best := -MaxDouble;
    for i := 0 to TestCount - 1 do
      if results[i, j] > best then
      begin
        best := results[i, j];
        frames[j].BandRMSWeighting := bwr[i];
        frames[j].BandBWeighting := not bwa[i];
        frames[j].BandWeightingApplyPower := bwp[i];
      end;
    writeln('#', j, #9, frames[j].BandRMSWeighting, #9, frames[j].BandBWeighting, #9, FloatToStr(frames[j].BandWeightingApplyPower));
  end;

  MakeFrames;
  MakeDstData;
end;

procedure TEncoder.SearchBestCorrRatios;

  function DoOne(ACR1, ACR2: Double; Verbose: Boolean): Double;
  begin
    CR1 := CICCorrRatio[BandDealiasSecondOrder, 0];
    CR2 := ACR1;

    PrepareFrames;
    MakeFrames;
    MakeDstData;
    SaveWAV;
    Result := -DoExternalEAQUAL(ParamStr(1), ParamStr(2), False, False, 2048);

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
{$if false}
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
      WriteLn('OutputBitDepth = ', enc.OutputBitDepth);
      WriteLn('MinChunkSize = ', enc.MinChunkSize);
      WriteLn('UsePython = ', BoolToStr(enc.UsePython, True));
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
      //enc.SaveRSC;

      enc.ComputeEAQUAL(enc.SampleCount, True, True, enc.srcData, enc.dstData);

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

