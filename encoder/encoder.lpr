program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, yakmo, ap, fft, conv, anysort;

const
  BandCount = 4;
  BandTransFactor = 0.1;
  LowCut = 30.0;
  HighCut = 16000.0;
  BandDealiasSecondOrder = True; // otherwise first order
  BandBWeighting = True; // otherwise A-weighting

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
    CorrectionFactor: Double;
    SecondOrder: Boolean;

    constructor Create(ARatio: Integer; ASecondOrder: Boolean);

    function ProcessSample(Smp: Double): Double;
  end;

  { TChunk }

  TChunk = class
  public
    encoder: TEncoder;
    band: TBand;
    reducedChunk: TChunk;

    index: Integer;

    mixList: TComplex1DArrayList;

    srcData: TDoubleDynArray;
    fft: TComplex1DArray;

    constructor Create(enc: TEncoder; bnd: TBand; idx: Integer);

    procedure ComputeFFT;
    procedure InitMix;
    procedure AddToMix(ch: TChunk);
    procedure FinalizeMix;
  end;

  TChunkList = specialize TFPGObjectList<TChunk>;

  { TBand }

  TBand = class
  public
    encoder: TEncoder;

    underSample, underSampleUnMin: Integer;
    chunkSize, chunkSizeUnMin: Integer;
    chunkCount: Integer;
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
    quality: Double;
    iterationCount: Integer;
    sampleRate: Integer;
    minChunkSize: Integer;
    srcDataCount: Integer;
    projectedDataCount: Integer;

    srcHeader: array[$00..$2b] of Byte;
    srcData: TDoubleDynArray;
    dstData: TSmallIntDynArray;

    bands: array[0..BandCount - 1] of TBand;

    class function make16BitSample(smp: Double): SmallInt;
    class function makeFloatSample(smp: SmallInt): Double;
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

    function DoFilterCoeffs(fc, transFactor: Double; HighPass: Boolean): TDoubleDynArray;
    function DoFilter(const samples, coeffs: TDoubleDynArray): TDoubleDynArray;
    function DoBPFilter(fcl, fch, transFactor: Double; underSample: Integer; const samples: TDoubleDynArray): TDoubleDynArray;

    function ComputeEAQUAL(chunkSz: Integer; UseDIX: Boolean; const smpRef, smpTst: TDoubleDynArray): Double;
    function ComputeEAQUALMulti(chunkSz: Integer; UseDIX: Boolean; const smpRef: TDoubleDynArray;
      smpTst: TDoubleDynArray2): TDoubleDynArray;
  end;


function IsDebuggerPresent () : LongBool stdcall; external 'kernel32.dll';

function CompareDouble(const d1, d2): integer;
var
  v1 : Double absolute d1;
  v2 : Double absolute d2;
begin
  Result := CompareValue(v1, v2);
end;

function GetToken(var s : string; const c : string) : string;
var
  p: Integer;
begin
  p := Pos(c, s);
  if p = 0 then p := Length(s);
  Result := Copy(s, 1, p-1);
  Delete(s, 1, p+Length(c)-1);
end;



constructor TSecondOrderCIC.Create(ARatio: Integer; ASecondOrder: Boolean);
begin
  SecondOrder := ASecondOrder;
  CorrectionFactor := ifthen(SecondOrder, 3.0, 1.0);
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

  hc := min(HighCut, encoder.sampleRate / 2);
  ratioP := log2(LowCut / hc) / BandCount + 0.5;
  ratioL := hc / encoder.sampleRate;

  if index = 0 then
    fcl := LowCut / encoder.sampleRate
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
  WriteLn('save #', index, ' ', fn);

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

  procedure DoChunk(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    chunk: TChunk;
  begin
    chunk := chunkList[AIndex];

    chunk.ComputeFFT;
  end;

var
  i: Integer;
  chunk: TChunk;
  cic: TSecondOrderCIC;
begin
  WriteLn('MakeChunks #', index, ' (', round(fcl * encoder.sampleRate), ' Hz .. ', round(fch * encoder.sampleRate), ' Hz); ', chunkSize, ' * (', chunkCount, ' -> ', desiredChunkCount,'); ', underSample);

  srcData := Copy(encoder.srcData);

  // compensate for decoder altering the pass band
  cic := TSecondOrderCIC.Create(underSample * 2, BandDealiasSecondOrder);
  try
    for i := -cic.Ratio + 1 to High(srcData) do
      srcData[Max(0, i)] := encoder.srcData[Max(0, i)] * (cic.CorrectionFactor + 1) -
        cic.CorrectionFactor * cic.ProcessSample(srcData[Min(High(srcData), i + cic.Ratio - 1)]);
  finally
    cic.Free;
  end;

  srcData := encoder.DoBPFilter(fcl, fch, BandTransFactor, 1, srcData);

  chunkList.Capacity := chunkCount;
  for i := 0 to chunkCount - 1 do
  begin
    chunk := TChunk.Create(encoder, Self, i);
    chunkList.Add(chunk);
  end;

  ProcThreadPool.DoParallelLocalProc(@DoChunk, 0, chunkCount - 1, nil);
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
  end;

var
  FN, Line: String;
  v1, v2, continuityFactor, sl, sh, fdl, fdh: Double;
  Dataset: TStringList;
  i, j, offset, itc: Integer;
  chunk: TChunk;
begin
  itc := encoder.iterationCount;
  if itc = 0 then Exit;
  if itc < 0 then itc := min(250, -itc * encoder.bands[BandCount - 1].chunkCount div chunkCount);

  WriteLn('KMeansReduce #', index, ' ', desiredChunkCount);

  FN := GetTempFileName('', 'dataset-'+IntToStr(GetCurrentThreadId)+'.txt');
  Dataset := TStringList.Create;
  Dataset.LineBreak := #10;

  continuityFactor := chunkSize / chunkSizeUnMin;
  offset := 4;
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

      Line := Line + Format('0:%e 1:%e 2:%e 3:%e ', [sl * continuityFactor, sh * continuityFactor, fdl * continuityFactor, fdh * continuityFactor]);

      for j := 0 to chunkSizeUnMin div 2 - 1 do
      begin
        v1 := chunkList[i].fft[j].X;
        v2 := chunkList[i].fft[j].Y;
        Line := Line + Format('%d:%e %d:%e ', [j * 2 + offset, v1, j * 2 + 1 + offset, v2]);
      end;
      Dataset.Add(Line);
    end;
    Dataset.SaveToFile(FN);
  finally
    Dataset.Free;
  end;

  SetLength(XYC, chunkList.Count);
  FillChar(XYC[0], chunkList.Count * SizeOF(Integer), $ff);
  DoExternalKMeans(FN, '', desiredChunkCount, itc, False, XYC);

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
  WriteLn('MakeDstData #', index);

  SetLength(dstData, Length(srcData));
  FillQWord(dstData[0], Length(srcData), 0);

  cic := TSecondOrderCIC.Create(underSample * 2, BandDealiasSecondOrder);
  pos := -cic.Ratio + 1;
  try
    for i := 0 to chunkList.Count - 1 do
    begin
      chunk := chunkList[i];

      for j := 0 to chunkSize - 1 do
      begin
        smp := chunk.reducedChunk.srcData[j];
        for k := 0 to underSample - 1 do
        begin
          dstData[max(0, pos)] := cic.ProcessSample(Smp);
          pos := min(High(dstData), pos + 1);
        end;
      end;
    end;

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

{ TEncoder }

procedure TEncoder.Load;
var
  fs: TFileStream;
  i: Integer;
begin
  WriteLn('load ', inputFN);
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

    bnd.KMeansReduce;
    bnd.MakeDstData;
    bnd.Save(ChangeFileExt(outputFN, '-' + IntToStr(AIndex) + '.wav'));
  end;

var
  i: Integer;
begin
  for i := 0 to BandCount - 1 do
    bands[i] := TBand.Create(Self, i);

  FindBestDesiredChunksCounts;

  for i := 0 to BandCount - 1 do
    bands[i].MakeChunks;

  ProcThreadPool.DoParallelLocalProc(@DoBand, 0, BandCount - 1, nil);
end;


procedure TEncoder.Save;
var
  fs: TFileStream;
begin
  WriteLn('save ', outputFN);

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
  fsq: Double;
  wgtCurve: array[0 .. BandCount - 1] of Double;
  i, allSz: Integer;
  sz: Double;
  full: Boolean;
begin
  for i := 0 to BandCount - 1 do
  begin
    bnd := bands[i];

    fsq := bnd.fcl * bnd.fch * sampleRate * sampleRate;

{$if BandBWeighting}
    // B-weighting
    wgtCurve[i] := sqr(12194.0) * sqrt(fsq) * fsq / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)) * sqrt(fsq + sqr(158.5)));
{$else}
    // A-weighting
    wgtCurve[i] := sqr(12194.0) * sqr(fsq) / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)) * sqrt((fsq + sqr(107.7)) * (fsq + sqr(737.9))));
{$endif}
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
  until (allSz >= srcDataCount * quality) or full;

  projectedDataCount := allSz;

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

function TEncoder.DoBPFilter(fcl, fch, transFactor: Double; underSample: Integer; const samples: TDoubleDynArray
  ): TDoubleDynArray;
var
  l, h: TDoubleDynArray;
begin
  l := DoFilterCoeffs(fcl * underSample, transFactor * (exp(fcl) - 1.0), True);
  h := DoFilterCoeffs(fch * underSample, transFactor * (exp(fch) - 1.0), False);

  Result := DoFilter(samples, l);
  Result := DoFilter(Result, h);
end;

constructor TEncoder.Create(InFN, OutFN: String);
begin
  inputFN := InFN;
  outputFN := OutFN;
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
  WriteLn('MakeDstData');

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

class function TEncoder.makeFloatSample(smp: SmallInt): Double;
begin
  Result := smp / -Low(SmallInt);
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
      WriteLn('Usage: (source file must be 16bit mono WAV)');
      writeln(ExtractFileName(ParamStr(0)) + ' <source file> <dest file> [quality 0.0-1.0] [iter count (0 = "lossless"; neg = relative to band)] [min chunk size 2-inf]');
      WriteLn;
      Exit;
    end;

    enc := TEncoder.Create(ParamStr(1), ParamStr(2));

    enc.quality := EnsureRange(StrToFloatDef(ParamStr(3), 0.5), 0.001, 1.0);
    enc.iterationCount := StrToIntDef(ParamStr(4), -10);
    enc.minChunkSize := StrToIntDef(ParamStr(5), 2);

    try

      enc.Load;
      enc.MakeBands;
      enc.MakeDstData;
      enc.Save;

    finally
      enc.Free;
    end;

{$if false}
    ShellExecute(0, 'open', PAnsiChar(ParamStr(2)), nil, nil, 0);
{$else}
    DoExternalEAQUAL(ParamStr(1), ParamStr(2), True, False, 2048);
{$endif}

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

