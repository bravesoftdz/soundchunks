program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, yakmo, ap, fft, conv, anysort, minlbfgs;

const
  BandCount = 4;
  CICCorrRatio: array[Boolean{2nd order?}, 0..1] of Double = (
    (0.51651, 0.47541),
    (2.05853, 2.26151)
  );

type
  TComplex1DArrayList = specialize TFPGList<TComplex1DArray>;
  TDoubleDynArray2 = array of TDoubleDynArray;

  TEncoder = class;
  TBand = class;

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
    encoder: TEncoder;
    band: TBand;
    reducedChunk, sameChunk: TChunk;

    index: Integer;
    bitRange: Integer;
    dstBitShift: Integer;
    tmpIndex: Integer;

    mixList: TComplex1DArrayList;

    srcData: TDoubleDynArray;
    fft: TComplex1DArray;
    dstData: TByteDynArray;

    constructor Create(enc: TEncoder; bnd: TBand; idx: Integer);

    procedure ComputeFFT;
    procedure InitMix;
    procedure AddToMix(ch: TChunk);
    procedure FinalizeMix;
    procedure ComputeBitRange;
    procedure MakeDstData;
  end;

  TChunkList = specialize TFPGObjectList<TChunk>;

  { TBand }

  TBand = class
  public
    encoder: TEncoder;

    underSample, underSampleUnMin: Integer;
    chunkSize, chunkSizeUnMin: Integer;
    chunkCount, uniqueChunkCount: Integer;
    index: Integer;
    fcl, fch: Double;

    desiredChunkCount: Integer;
    srcData: TDoubleDynArray;
    dstData: TDoubleDynArray;

    chunkList: TChunkList;
    reducedChunks: TChunkList;

    constructor Create(enc: TEncoder; idx: Integer);
    destructor Destroy; override;

    procedure Save(fn: String);

    procedure MakeChunks;
    procedure KMeansReduce;
    procedure MakeDstData;
  end;

  { TEncoder }

  TEncoder = class
  public
    inputFN, outputFN: String;

    Quality: Double;
    Precision: Integer;
    BandTransFactor: Double;
    LowCut: Double;
    HighCut: Double;
    BandDealiasSecondOrder: Boolean; // otherwise first order
    BandBWeighting: Boolean; // otherwise A-weighting
    BandWeightingApplyCount: Integer;
    OutputBitDepth: Integer; // max 8Bits
    UseCUDA: Boolean;

    sampleRate: Integer;
    minChunkSize: Integer;
    srcDataCount: Integer;
    projectedDataCount: Integer;
    verbose: Boolean;

    CRSearch: Boolean;
    CR1, CR2: Double;

    srcHeader: array[$00..$2b] of Byte;
    srcData: TDoubleDynArray;
    dstData: TSmallIntDynArray;

    bands: array[0..BandCount - 1] of TBand;

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
    procedure Save;

    procedure FindBestDesiredChunksCounts;
    procedure MakeBands;
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


constructor TBand.Create(enc: TEncoder; idx: Integer);
var
  hc, ratioP, ratioL: Double;
begin
  encoder := enc;
  index := idx;

  // determing low and high bandpass frequencies

  hc := min(encoder.HighCut, encoder.sampleRate / 2);
  ratioP := log2(encoder.LowCut / hc) / BandCount + 0.5;
  ratioL := hc / encoder.sampleRate;

  if index = 0 then
    fcl := encoder.LowCut / encoder.sampleRate
  else
    fcl := power(2.0, (BandCount - index) * ratioP) * ratioL;

  if index = BandCount - 1 then
    fch := hc / encoder.sampleRate
  else
    fch := power(2.0, (BandCount - 1 - index) * ratioP) * ratioL;

  // undersample if the band high freq is a lot lower than nyquist

  chunkSize := round(intpower(2.0, floor(-log2(fcl))));
  underSample := round(intpower(2.0, floor(-log2(fch)) - 2));
  underSample := Max(1, underSample);
  chunkSize := chunkSize div underSample;

  underSampleUnMin := underSample;
  chunkSizeUnMin := chunkSize;
  if chunkSize < encoder.minChunkSize then
  begin
    underSample := max(1, (underSample * chunkSize) div encoder.minChunkSize);
    chunkSize := encoder.minChunkSize;
  end;

  chunkCount := (encoder.srcDataCount - 1) div (chunkSize * underSample) + 1;
  uniqueChunkCount := chunkCount;

  chunkList := TChunkList.Create;
  reducedChunks := TChunkList.Create;
end;

destructor TBand.Destroy;
begin
  chunkList.Free;
  reducedChunks.Free;

  inherited Destroy;
end;

procedure TBand.Save(fn: String);
var
  i: Integer;
  fs: TFileStream;
begin
  WriteLn('Save #', index, ' ', fn);

  fs := TFileStream.Create(fn, fmCreate or fmShareDenyWrite);
  try
    fs.WriteBuffer(encoder.srcHeader[0], SizeOf(encoder.srcHeader));
    for i := 0 to encoder.srcDataCount - 1 do
      fs.WriteWord(Word(TEncoder.make16BitSample(dstData[i])));
  finally
    fs.Free;
  end;
end;

procedure TBand.MakeChunks;
var
  i: Integer;
  chunk: TChunk;
  cic: TSecondOrderCIC;
  cicSmp, cr1, cr2: Double;
begin
  if not encoder.CRSearch then
    WriteLn('MakeChunks #', index, ' (', round(fcl * encoder.sampleRate), ' Hz .. ', round(fch * encoder.sampleRate), ' Hz); ', chunkSize, ' * (', chunkCount, ' -> ', desiredChunkCount,'); ', underSample);

  srcData := encoder.DoBPFilter(fcl, fch, encoder.BandTransFactor, encoder.srcData);

  if underSample > 1 then
  begin
    // compensate for decoder altering the pass band
    cic := TSecondOrderCIC.Create(underSample * 2, encoder.BandDealiasSecondOrder);

    cr1 := CICCorrRatio[encoder.BandDealiasSecondOrder, 0];
    cr2 := CICCorrRatio[encoder.BandDealiasSecondOrder, 1];
    if encoder.CRSearch then
    begin
      cr1 := encoder.CR1;
      cr2 := encoder.CR2;
    end;

    try
      for i := -cic.Ratio + 1 to High(srcData) do
      begin
        cicSmp := cic.ProcessSample(srcData[Min(High(srcData), i + cic.Ratio - 1)]);
        if i >= 0 then
          srcData[i] += srcData[i] * cr1 - cicSmp * cr2;
      end;
    finally
      cic.Free;
    end;
  end;

  chunkList.Capacity := chunkCount;
  for i := 0 to chunkCount - 1 do
  begin
    chunk := TChunk.Create(encoder, Self, i);
    if not encoder.CRSearch then
      chunk.ComputeFFT;
    chunk.MakeDstData;
    chunkList.Add(chunk);
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

    for i := 0 to chunkList.Count - 1 do
      if XYC[i] = AIndex then
      begin
        reducedChunk.AddToMix(chunkList[i]);
        chunkList[i].reducedChunk := reducedChunk;
      end;

    reducedChunk.FinalizeMix;
    reducedChunk.MakeDstData;
  end;

var
  FN, Line: String;
  v1, v2, continuityFactor, sl, sh, fdl, fdh: Double;
  Dataset: TStringList;
  i, j, prec: Integer;
  chunk: TChunk;
begin
  prec := encoder.Precision;
  if (prec = 0) or (desiredChunkCount = chunkCount) then Exit;

  WriteLn('KMeansReduce #', index, ' ', desiredChunkCount);

  FN := GetTempFileName('', 'dataset-'+IntToStr(GetCurrentThreadId)+'.txt');
  Dataset := TStringList.Create;
  Dataset.LineBreak := #10;

  continuityFactor := chunkSize / chunkSizeUnMin;
  try
    for i := 0 to chunkList.Count - 1 do
    begin
      Line := IntToStr(i) + ' ';

      // add first and last sample to features to allow continuity between chunks
      sl := chunkList[i].srcData[0];
      sh := chunkList[i].srcData[chunkSize - 1];
      // add approximate derivatives of first and last sample to features to allow continuity between chunks
      fdl :=  -chunkList[i].srcData[0] + chunkList[i].srcData[1];
      fdh :=  chunkList[i].srcData[chunkSize - 1] - chunkList[i].srcData[chunkSize - 2];

      Line := Line + Format('%0.8f %0.8f %0.8f %0.8f ', [sl * continuityFactor, sh * continuityFactor, fdl * continuityFactor, fdh * continuityFactor]);

      for j := 0 to chunkSizeUnMin - 1 do
      begin
        v1 := chunkList[i].fft[j].X;
        v2 := chunkList[i].fft[j].Y;
        Line := Line + Format('%0.8f %0.8f ', [v1, v2]);
      end;
      Dataset.Add(Line);
    end;
    Dataset.SaveToFile(FN);
  finally
    Dataset.Free;
  end;

  SetLength(XYC, chunkList.Count);
  FillChar(XYC[0], chunkList.Count * SizeOF(Integer), $ff);
  DoExternalKMeans(FN, desiredChunkCount, 1, prec, encoder.verbose, encoder.UseCUDA, XYC);

  for i := 0 to desiredChunkCount - 1 do
  begin
    chunk := TChunk.Create(encoder, Self, 0);

    reducedChunks.Add(chunk);
    DoXYC(i);
  end;
end;

procedure TBand.MakeDstData;
var
  i, j, k, pos: Integer;
  chunk: TChunk;
  smp: Double;
  cic: TSecondOrderCIC;
begin
  //WriteLn('MakeDstData #', index);

  SetLength(dstData, Length(srcData));
  FillQWord(dstData[0], Length(srcData), 0);

  if underSample > 1 then
  begin
    cic := TSecondOrderCIC.Create(underSample * 2, encoder.BandDealiasSecondOrder);
    pos := -cic.Ratio + 1;
  end
  else
  begin
    cic := nil;
    pos := 0;
  end;

  try
    for i := 0 to chunkList.Count - 1 do
    begin
      chunk := chunkList[i];

      for j := 0 to chunkSize - 1 do
      begin
        smp := TEncoder.makeFloatSample(chunk.reducedChunk.dstData[j], encoder.OutputBitDepth, chunk.reducedChunk.dstBitShift);
        for k := 0 to underSample - 1 do
        begin
          if Assigned(cic) then
            dstData[max(0, pos)] := cic.ProcessSample(Smp)
          else
            dstData[pos] := Smp;

          pos := min(High(dstData), pos + 1);
        end;
      end;
    end;

    if Assigned(cic) then
      for k := 0 to underSample - 1 do
      begin
        dstData[pos] := cic.ProcessSample(0);
        pos := min(High(dstData), pos + 1);
      end;
  finally
    cic.Free;
  end;
end;

{ TChunk }

constructor TChunk.Create(enc: TEncoder; bnd: TBand; idx: Integer);
var
  j, k, pos: Integer;
  acc: Double;
begin
  index := idx;
  encoder := enc;
  band := bnd;

  SetLength(srcData, band.chunkSize);

  for j := 0 to band.chunkSize - 1 do
  begin
    pos := (idx * band.chunkSize + j) * band.underSample;
    acc := 0.0;
    for k := 0 to band.underSample - 1 do
      acc += band.srcData[EnsureRange(pos - k, 0, High(band.srcData))];
    srcData[j] := acc / band.underSample;
  end;

  reducedChunk := Self;
end;

procedure TChunk.ComputeFFT;
begin
  FFTR1D(srcData, band.chunkSize, fft);
end;

procedure TChunk.InitMix;
begin
  mixList := TComplex1DArrayList.Create;
end;

procedure TChunk.AddToMix(ch: TChunk);
begin
  mixList.Add(ch.fft);
end;

procedure TChunk.FinalizeMix;
var
  i, k: Integer;
  acc: Complex;
begin
  if not Assigned(mixList) or (mixList.Count = 0) then
    Exit;

  SetLength(fft, band.chunkSize);
  for k := 0 to band.chunkSize - 1 do
  begin
    acc.X := 0.0;
    acc.Y := 0.0;

    for i := 0 to mixList.Count - 1 do
    begin
      acc.X += mixList[i][k].X;
      acc.Y += mixList[i][k].Y;
    end;

    acc.X /= mixList.Count;
    acc.Y /= mixList.Count;

    fft[k] := acc;
  end;

  FFTR1DInv(fft, band.chunkSize, srcData);

  FreeAndNil(mixList);
end;

procedure TChunk.ComputeBitRange;
var
  i, hiSmp: Integer;
begin
  hiSmp := 0;
  for i := 0 to band.chunkSize - 1 do
    hiSmp := max(hiSmp, abs(TEncoder.make16BitSample(srcData[i])));
  bitRange := EnsureRange(1 + ceil(log2(hiSmp + 1.0)), encoder.OutputBitDepth, 16);
end;

procedure TChunk.MakeDstData;
var
  i: Integer;
begin
  ComputeBitRange;
  dstBitShift := bitRange - encoder.OutputBitDepth;
  SetLength(dstData, band.chunkSize);
  for i := 0 to band.chunkSize - 1 do
    dstData[i] := TEncoder.makeOutputSample(srcData[i], encoder.OutputBitDepth, dstBitShift);
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
    srcDataCount := (fs.Size - fs.Position) div 2;
    SetLength(srcData, srcDataCount);
    FillQWord(srcData[0], srcDataCount, 0);
    for i := 0 to srcDataCount - 1 do
      srcData[i] := makeFloatSample(SmallInt(fs.ReadWord));
  finally
    fs.Free;
  end;

  sampleRate := PInteger(@srcHeader[$18])^;
  writeln('sampleRate = ', sampleRate);
  WriteLn('srcDataCount = ', srcDataCount);
end;

procedure TEncoder.MakeBands;

  procedure DoBand(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    bnd: TBand;
  begin
    bnd := bands[AIndex];

    if not CRSearch then
      bnd.KMeansReduce;

    bnd.MakeDstData;

    if not CRSearch then
      bnd.Save(ChangeFileExt(outputFN, '-' + IntToStr(AIndex) + '.wav'));
  end;

  procedure DoBandChunks(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    bnd: TBand;
  begin
    bnd := bands[AIndex];

    bnd.MakeChunks;
  end;

var
  i: Integer;
begin
  for i := 0 to BandCount - 1 do
    bands[i] := TBand.Create(Self, i);

  FindBestDesiredChunksCounts;

  ProcThreadPool.DoParallelLocalProc(@DoBandChunks, 0, BandCount - 1, nil);

  if UseCUDA then
  begin
    for i := 0 to BandCount - 1 do
      DoBand(i, nil, nil);
  end
  else
  begin
    ProcThreadPool.DoParallelLocalProc(@DoBand, 0, BandCount - 1, nil);
  end;
end;


procedure TEncoder.Save;
var
  fs: TFileStream;
begin
  if not CRSearch then
    WriteLn('Save ', outputFN);

  fs := TFileStream.Create(outputFN, fmCreate or fmShareDenyWrite);
  try
    fs.WriteBuffer(srcHeader[0], SizeOf(srcHeader));
    fs.WriteBuffer(dstData[0], srcDataCount * 2);
  finally
    fs.Free;
  end;
end;

procedure TEncoder.FindBestDesiredChunksCounts;
var
  bnd: TBand;
  wgt, fsq: Double;
  wgtCurve: array[0 .. BandCount - 1] of Double;
  i, allSz: Integer;
  sz: Double;
  full: Boolean;
begin
  for i := 0 to BandCount - 1 do
  begin
    bnd := bands[i];

    fsq := bnd.fcl * bnd.fch * sampleRate * sampleRate;

    if BandBWeighting then
      // B-weighting
      wgt := sqr(12194.0) * sqrt(fsq) * fsq / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)) * sqrt(fsq + sqr(158.5)))
    else
      // A-weighting
      wgt := sqr(12194.0) * sqr(fsq) / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)) * sqrt((fsq + sqr(107.7)) * (fsq + sqr(737.9))));

    wgtCurve[i] := intpower(wgt, BandWeightingApplyCount);
  end;

  sz := 1;
  repeat
    full := True;
    allSz := 0;
    for i := 0 to BandCount - 1 do
    begin
      bnd := bands[i];

      bnd.desiredChunkCount := min(bnd.chunkCount, round(sz * wgtCurve[i]));
      full := full and (bnd.desiredChunkCount = bnd.chunkCount);

      allSz += bnd.desiredChunkCount * bnd.chunkSize;
    end;
    sz += 1.0;
  until (allSz >= srcDataCount * Quality) or full;

  projectedDataCount := allSz;

  if not CRSearch then
    WriteLn('projectedDataCount = ', projectedDataCount);
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

  Quality := 0.5;
  Precision := 5;
  BandTransFactor := 0.05;
  LowCut := 30.0;
  HighCut := 18000.0;
  BandDealiasSecondOrder := True;
  BandBWeighting := True;
  BandWeightingApplyCount := 1;
  OutputBitDepth := 8;
  UseCUDA := False;
end;

destructor TEncoder.Destroy;
var
  i: Integer;
begin
  for i := 0 to BandCount - 1 do
    bands[i].Free;

  inherited Destroy;
end;

procedure TEncoder.MakeDstData;
var
  i, j: Integer;
  acc: Double;
begin
  //WriteLn('MakeDstData');

  SetLength(dstData, Length(srcData));
  FillWord(dstData[0], Length(srcData), 0);
  for i := 0 to High(dstData) do
  begin
    acc := 0.0;
    for j := 0 to BandCount - 1 do
      acc += bands[j].dstData[i];
    dstData[i] := make16BitSample(acc);
  end;
end;

procedure TEncoder.SearchBestCorrRatios;

var
  srcf, dstf: TDoubleDynArray;

  function DoOne(ACR1, ACR2: Double; Verbose: Boolean): Double;
  var
    i: Integer;
    eaq: Double;
  begin
    CR1 := ACR1;
    CR2 := ACR2;

    MakeBands;
    MakeDstData;

    for i := 0 to srcDataCount - 1 do
      dstf[i] := makeFloatSample(dstData[i]);

    Result := CompareDCT(0, srcDataCount - 1, False, srcf, dstf);

    if Verbose then
    begin
      Save;
      eaq := DoExternalEAQUAL(ParamStr(1), ParamStr(2), False, False, 2048);
      Write(#13'(', FormatFloat('0.00000', CR1), ', ', FormatFloat('0.00000', CR2), ')'#9, FormatFloat('0.00000', Result),#9, FloatToStr(eaq), '                ');
    end;
  end;

var
  N : AlglibInteger;
  State : MinLBFGSState;
  Rep : MinLBFGSReport;
  S : TReal1DArray;
  H : Double;
begin
  srcf := DoBPFilter(bands[0].fcl, bands[BandCount - 1].fch, BandTransFactor, srcData);
  SetLength(dstf, srcDataCount);

  CRSearch := True;

  H := 1e-3;
  N := 2;
  SetLength(S, 2);
  S[0] := 1.0 + 2.0 * Ord(BandDealiasSecondOrder);
  S[1] := S[0];
  MinLBFGSCreate(N, N, S, State);
  MinLBFGSSetCond(State, H, 0.0, 0.0, 0);
  MinLBFGSSetStpMax(State, 0.1);
  while MinLBFGSIteration(State) do
  begin
    if State.NeedFG then
    begin
      State.F := DoOne(State.X[0], State.X[1], True);

      State.G[0] := (
          -DoOne(State.X[0] + 2 * H, State.X[1], False) +
          8 * DoOne(State.X[0] + H, State.X[1], False) +
          -8 * DoOne(State.X[0] - H, State.X[1], False) +
          DoOne(State.X[0] - 2 * H, State.X[1], False)
        ) / (12 * H);

      State.G[1] := (
          -DoOne(State.X[0], State.X[1] + 2 * H, False) +
          8 * DoOne(State.X[0], State.X[1] + H, False) +
          -8 * DoOne(State.X[0], State.X[1] - H, False) +
          DoOne(State.X[0], State.X[1] - 2 * H, False)
        ) / (12 * H);
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

begin
  try
    FormatSettings.DecimalSeparator := '.';

{$ifdef DEBUG}
    ProcThreadPool.MaxThreadCount := 1;
{$endif}

    if ParamCount < 2 then
    begin
      WriteLn('Usage: ', ExtractFileName(ParamStr(0)) + ' <source file> <dest file> [options]');
      WriteLn(#9'-q'#9'encoder quality (0.0-1.0); example: "-q0.2"');
      WriteLn(#9'-pr'#9'K-means precision (0-9); 0: "lossless" mode');
      WriteLn(#9'-btf'#9'band transition factor (0.0001-1)');
      WriteLn(#9'-lc'#9'bass cutoff frequency');
      WriteLn(#9'-hc'#9'treble cutoff frequency');
      WriteLn(#9'-bd1'#9'use first order dealias filter (otherwise second order)');
      WriteLn(#9'-bwa'#9'use A-weighting for bands (otherwise B-weighting)');
      WriteLn(#9'-bwc'#9'band weighting apply count (0-inf)');
      WriteLn(#9'-obd'#9'output bit depth (1-8)');
      WriteLn(#9'-v'#9'verbose K-means');
      WriteLn(#9'-cu'#9'use CUDA GPU K-means');
      WriteLn;
      Writeln('(source file must be 16bit mono WAV)');
      WriteLn;
      Exit;
    end;

    enc := TEncoder.Create(ParamStr(1), ParamStr(2));
    try
      enc.Quality :=  ParamValue('-q', enc.Quality);
      enc.Precision := round(ParamValue('-pr', enc.Precision));
      enc.BandTransFactor :=  ParamValue('-btf', enc.BandTransFactor);
      enc.LowCut :=  ParamValue('-lc', enc.LowCut);
      enc.HighCut :=  ParamValue('-hc', enc.HighCut);
      enc.BandDealiasSecondOrder :=  ParamStart('-bd1') = -1;
      enc.BandBWeighting :=  ParamStart('-bwa') = -1;
      enc.BandWeightingApplyCount := round(ParamValue('-bwc', enc.BandWeightingApplyCount));
      enc.OutputBitDepth :=  round(ParamValue('-obd', enc.OutputBitDepth));
      enc.verbose := ParamStart('-v') <> -1;
      enc.UseCUDA := ParamStart('-cu') <> -1;

      WriteLn('Quality = ', FloatToStr(enc.Quality));
      WriteLn('Precision = ', enc.Precision);
      WriteLn('BandTransFactor = ', FloatToStr(enc.BandTransFactor));
      WriteLn('LowCut = ', FloatToStr(enc.LowCut));
      WriteLn('HighCut = ', FloatToStr(enc.HighCut));
      WriteLn('BandDealiasSecondOrder = ', BoolToStr(enc.BandDealiasSecondOrder, True));
      WriteLn('BandBWeighting = ', BoolToStr(enc.BandBWeighting, True));
      WriteLn('BandWeightingApplyCount = ', enc.BandWeightingApplyCount);
      WriteLn('OutputBitDepth = ', enc.OutputBitDepth);
      WriteLn('UseCUDA = ', BoolToStr(enc.UseCUDA, True));
      WriteLn;

      enc.Load;
      enc.MakeBands;
      enc.MakeDstData;
      enc.Save;

      DoExternalEAQUAL(ParamStr(1), ParamStr(2), True, True, 2048);

      if ParamStart('-scr') <> -1 then
        enc.SearchBestCorrRatios;

    finally
      enc.Free;
    end;

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

