program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, yakmo, ap, conv;

const
  BandCount = 4;
  BandTransFactor = 0.1;
  LowCut = 30.0;
  HighCut = 18000.0;

  Butter2ndCoeffs: array[0..9, 0..2] of Double =
  (
    (0.0000000000, -0.1715728753, 3.414213562e+00),
    (0.9428090416, -0.3333333333, 1.024264069e+01),
    (1.4542435863, -0.5740619151, 3.338387406e+01),
    (1.7237761728, -0.7575469445, 1.184456202e+02),
    (1.8613611468, -0.8703674775, 4.441320406e+02),
    (1.9306064272, -0.9329347318, 1.717988320e+03),
    (1.9652933726, -0.9658854606, 6.755753138e+03),
    (1.9826454185, -0.9827947193, 2.679155178e+04),
    (1.9913225484, -0.9913600352, 1.067042553e+05),
    (1.9956612539, -0.9956706459, 4.258941024e+05)
  );

type
  TDoubleDynArray2 = array of TDoubleDynArray;

  TEncoder = class;
  TBand = class;

  { TChunk }

  TChunk = class
  public
    encoder: TEncoder;
    band: TBand;
    reducedChunk: TChunk;

    index: Integer;
    dctAddCount: Integer;

    srcData: TDoubleDynArray;
    dct: TDoubleDynArray;

    constructor Create(enc: TEncoder; bnd: TBand; idx: Integer);

    procedure ComputeDCT;
    procedure InitDCTAdd;
    procedure AddToDCT(ch: TChunk);
    procedure FinalizeDCTAdd;
  end;

  TChunkList = specialize TFPGObjectList<TChunk>;

  { TBand }

  TBand = class
  public
    encoder: TEncoder;

    underSample: Integer;
    chunkSize: Integer;
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
    quality: Double;
    restartCount: Integer;
    sampleRate: Integer;
    minChunkSize: Integer;
    srcDataCount: Integer;
    projectedDataCount: Integer;

    srcHeader: array[$00..$2b] of Byte;
    srcData: TDoubleDynArray;
    dstData: TSmallIntDynArray;

    bands: array[0..BandCount - 1] of TBand;

    class function make16BitSample(smp: Double): SmallInt;
    class function ComputeDCT(chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
    class function ComputeInvDCT(chunkSz: Integer; const dct: TDoubleDynArray): TDoubleDynArray;
    class function CompareDCT(firstCoeff, lastCoeff: Integer; compress: Boolean; const dctA, dctB: TDoubleDynArray): Double;
    class function CompressDCT(coeff: Double): Double;
    class function CheckJoinPenalty(x, y, z, a, b, c: Double; TestRange: Boolean): Boolean; inline;
    class function DoButterWorthLP(fcDiv: Integer; x, xm1, xm2, ym1, ym2: Double): Double;

    constructor Create;
    destructor Destroy; override;

    procedure Load(fn: String);
    procedure Save(fn: String);

    procedure FindBestDesiredChunksCounts;
    procedure MakeBands;
    procedure MakeDstData;

    function DoFilterCoeffs(fc, transFactor: Double; HighPass: Boolean): TDoubleDynArray;
    function DoFilter(fc, transFactor: Double; HighPass: Boolean; const samples: TDoubleDynArray): TDoubleDynArray;
    function DoBPFilter(fcl, fch, transFactor: Double; chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;

    function ComputeEAQUAL(chunkSz: Integer; UseDIX: Boolean; const smpRef, smpTst: TDoubleDynArray): Double;
    function ComputeEAQUALMulti(chunkSz: Integer; UseDIX: Boolean; const smpRef: TDoubleDynArray;
      smpTst: TDoubleDynArray2): TDoubleDynArray;
  end;


constructor TBand.Create(enc: TEncoder; idx: Integer);
var
  ratio: Double;
  fcdc: Double;
begin
  encoder := enc;
  index := idx;

  fcdc := LowCut / encoder.sampleRate;
  ratio := log2(HighCut / LowCut) / BandCount;
  fcl := fcdc * power(2.0, index * ratio);
  fch := fcdc * power(2.0, (index + 1) * ratio);

  chunkSize := round(intpower(2.0, ceil(-log2(fcl))));
  underSample := round(intpower(2.0, floor(-log2(fch))));
  chunkSize := chunkSize div underSample;

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

    chunk.ComputeDCT;
  end;

var
  i: Integer;
  chunk: TChunk;
begin
  WriteLn('MakeChunks #', index, ' (', round(fcl * encoder.sampleRate), ' Hz .. ', round(fch * encoder.sampleRate), ' Hz); ', chunkSize, ' * (', chunkCount, ' -> ', desiredChunkCount,')');

  srcData := encoder.DoBPFilter(fcl, fch, BandTransFactor, chunkSize, encoder.srcData);

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
    i, First: Integer;
    reducedChunk: TChunk;
  begin
    First := -1;
    for i := 0 to chunkList.Count - 1 do
      if XYC[i] = AIndex then
      begin
        First := i;
        Break;
      end;

    reducedChunk := reducedChunks[AIndex];
    reducedChunk.InitDCTAdd;

    if First <> -1 then
      for i := 0 to chunkList.Count - 1 do
        if XYC[i] = AIndex then
        begin
          reducedChunk.AddToDCT(chunkList[i]);
          chunkList[i].reducedChunk := reducedChunk;
        end;

    reducedChunk.FinalizeDCTAdd;
  end;

var
  FN, Line: String;
  v1: Double;
  Dataset: TStringList;
  i, j : Integer;
begin
  exit;

  WriteLn('KMeansReduce #', index, ' ', desiredChunkCount);

  FN := GetTempFileName('', 'dataset-'+IntToStr(GetCurrentThreadId)+'.txt');
  Dataset := TStringList.Create;
  Dataset.LineBreak := #10;

  try
    for i := 0 to chunkList.Count - 1 do
    begin
      Line := IntToStr(i) + ' ';
      for j := 0 to chunkSize - 1 do
      begin
        v1 := chunkList[i].dct[j];
        //v1 := TEncoder.CompressDCT(v1);
        Line := Line + Format('%d:%.12g ', [j, v1]);
      end;
      Dataset.Add(Line);
    end;
    Dataset.SaveToFile(FN);
  finally
    Dataset.Free;
  end;

  SetLength(XYC, chunkList.Count);
  FillChar(XYC[0], chunkList.Count * SizeOF(Integer), $ff);
  DoExternalKMeans(FN, desiredChunkCount, encoder.RestartCount, False, XYC);

  for i := 0 to desiredChunkCount - 1 do
  begin
    reducedChunks.Add(TChunk.Create(encoder, Self, 0));
    DoXYC(i);
  end;
end;

{$if true}
procedure TBand.MakeDstData;
var
  i, j, k: Integer;
  chunk: TChunk;
  v, vv, smp, alpha: Double;
begin
  WriteLn('MakeDstData #', index);

  v := 0;
  vv := 0;
  SetLength(dstData, Length(srcData));
  FillQWord(dstData[0], Length(srcData), 0);

  alpha := 2.0 / underSample;
  alpha := alpha / (alpha + 1);

  for i := 0 to chunkList.Count - 1 do
  begin
    chunk := chunkList[i];

    for j := 0 to chunkSize - 1 do
    begin
      smp := chunk.reducedChunk.srcData[j];

      for k := 0 to underSample - 1 do
      begin
        v := (1.0 - alpha) * v + alpha * smp;
        vv := (1.0 - alpha) * vv + alpha * v;

        dstData[(i * chunkSize + j) * underSample + k] := vv;
      end;
    end;
  end;
end;
{$else}
procedure TBand.MakeDstData;
var
  i, j, k: Integer;
  chunk: TChunk;
  smp, psmp, ppsmp, dst, pdst, ppdst: Double;
begin
  WriteLn('MakeDstData #', index);

  SetLength(dstData, Length(srcData));
  FillQWord(dstData[0], Length(srcData), 0);

  psmp := 0;
  ppsmp := 0;
  pdst := 0;
  ppdst := 0;
  for i := 0 to chunkList.Count - 1 do
  begin
    chunk := chunkList[i];

    for j := 0 to chunkSize - 1 do
    begin
      smp := chunk.reducedChunk.srcData[j];
      if underSample = 1 then
        dstData[i * chunkSize + j] := smp
      else
        for k := 0 to underSample - 1 do
        begin
          dst := TEncoder.DoButterWorthLP(underSample, smp, psmp, ppsmp, pdst, ppdst);

          dstData[(i * chunkSize + j) * underSample + k] := dst;

          ppdst := pdst;
          pdst := dst;

          ppsmp := psmp;
          psmp := smp;
        end;

    end;
  end;
end;
{$endif}

{ TChunk }

constructor TChunk.Create(enc: TEncoder; bnd: TBand; idx: Integer);
var
  j: Integer;
begin
  index := idx;
  encoder := enc;
  band := bnd;

  SetLength(srcData, band.chunkSize);

  for j := 0 to band.chunkSize - 1 do
    srcData[j] := band.srcData[(idx * band.chunkSize + j) * band.underSample];

  reducedChunk := Self;
  SetLength(dct, band.chunkSize);
end;

procedure TChunk.ComputeDCT;
begin
  dct := TEncoder.ComputeDCT(band.chunkSize, srcData);
end;

procedure TChunk.InitDCTAdd;
begin
  dctAddCount := 0;
  FillQWord(dct[0], band.chunkSize, 0);
end;

procedure TChunk.AddToDCT(ch: TChunk);
var
  k: Integer;
begin
  for k := 0 to band.chunkSize - 1 do
    dct[k] += ch.dct[k];
  Inc(dctAddCount);
end;

procedure TChunk.FinalizeDCTAdd;
var
  k: Integer;
begin
  if dctAddCount = 0 then Exit;

  for k := 0 to band.chunkSize - 1 do
    dct[k] /= dctAddCount;

  srcData := TEncoder.ComputeInvDCT(band.chunkSize, dct);
end;

{ TEncoder }

procedure TEncoder.Load(fn: String);
var
  fs: TFileStream;
  i: Integer;
begin
  WriteLn('load ', fn);
  fs := TFileStream.Create(fn, fmOpenRead or fmShareDenyNone);
  try
    fs.ReadBuffer(srcHeader[0], SizeOf(srcHeader));
    srcDataCount := (fs.Size - fs.Position) div 2;
    SetLength(srcData, srcDataCount + 65536);
    FillQWord(srcData[0], srcDataCount + 65536, 0);
    for i := 0 to srcDataCount - 1 do
      srcData[i] := SmallInt(fs.ReadWord);
  finally
    fs.Free;
  end;

  sampleRate := PInteger(@srcHeader[$18])^;
  writeln(sampleRate, ' Hz');
end;

procedure TEncoder.MakeBands;

  procedure DoBand(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    bnd: TBand;
  begin
    bnd := bands[AIndex];

    bnd.KMeansReduce;
    bnd.MakeDstData;
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


procedure TEncoder.Save(fn: String);
var
  i: Integer;
  fs: TFileStream;
begin
  WriteLn('save ', fn);

  for i := 0 to BandCount - 1 do
    bands[i].Save(ChangeFileExt(fn, '-' + IntToStr(i) + '.wav'));

  fs := TFileStream.Create(fn, fmCreate or fmShareDenyWrite);
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
  fsq, dbBatt: Double;
  i, sz, allSz: Integer;
begin
  sz := 1;
  repeat
    allSz := 0;
    for i := 0 to BandCount - 1 do
    begin
      bnd := bands[i];

      fsq := bnd.fcl * sampleRate * bnd.fch * sampleRate;
      dbBatt := sqr(12194.0) * fsq * sqrt(fsq) / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)) * sqrt(fsq + sqr(158.5)));

      bnd.desiredChunkCount := min(bnd.chunkCount, round(sz * dbBatt));
      allSz += bnd.desiredChunkCount * bnd.chunkSize;
    end;
    Inc(sz);
  until allSz >= srcDataCount * quality;

  projectedDataCount := allSz;

  WriteLn('projectedDataCount = ', projectedDataCount);
end;

function TEncoder.DoFilter(fc, transFactor: Double; HighPass: Boolean; const samples: TDoubleDynArray
  ): TDoubleDynArray;
var
  i: Integer;
  h: TReal1DArray;
begin
  if (fc <= 0.0) and HighPass or (fc >=0.5) and not HighPass then
  begin
    Result := samples;
    Exit;
  end;

  h := DoFilterCoeffs(fc, transFactor, HighPass);

  ConvR1D(samples, Length(samples), h, Length(h), Result);

  for i := 0 to High(samples) do
    Result[i] := Result[i + (Length(h) - 1) div 2];

  SetLength(Result, Length(samples));
end;

function TEncoder.DoBPFilter(fcl, fch, transFactor: Double; chunkSz: Integer; const samples: TDoubleDynArray
  ): TDoubleDynArray;
begin
  Assert(fch > fcl);
  Result := DoFilter(fcl, transFactor, True, samples);
  Result := DoFilter(fch, transFactor, False, Result);
end;

constructor TEncoder.Create;
begin

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

class function TEncoder.DoButterWorthLP(fcDiv: Integer; x, xm1, xm2, ym1, ym2: Double): Double;
begin
  fcDiv := BsrDWord(fcDiv);
  Result := (xm2 + 2.0 * xm1 + x) / Butter2ndCoeffs[fcDiv][2] + Butter2ndCoeffs[fcDiv][1] * ym2 + Butter2ndCoeffs[fcDiv][0] * ym1;
end;

function TEncoder.DoFilterCoeffs(fc, transFactor: Double; HighPass: Boolean): TDoubleDynArray;
var
  b, sinc, win, sum: Double;
  i, N: Integer;
begin
  b := fc * transFactor;
  N := ceil(4 / b);
  if (N mod 2) = 0 then N += 1;

  //writeln('DoFilter ', ifthen(HighPass, 'HP', 'LP'), ' ', FloatToStr(sampleRate * fc), ' ', N);

  SetLength(Result, N);
  sum := 0;
  for i := 0 to N - 1 do
  begin
    sinc := 2.0 * fc * (i - (N - 1) / 2.0) * pi;
    if IsZero(sinc) then
      sinc := 1.0
    else
      sinc := sin(sinc) / sinc;

    win := 0.42 - 0.5 * cos(2 * pi * i / (N - 1)) + 0.08 * cos(4 * pi * i / (N - 1));

    Result[i] := sinc * win;
    sum += Result[i];
  end;

  if HighPass then
  begin
    for i := 0 to N - 1 do
      Result[i] := -Result[i] / sum;

    Result[(N - 1) div 2] += 1;
  end
  else
  begin
    for i := 0 to N - 1 do
      Result[i] := Result[i] / sum;
  end;
end;

class function TEncoder.make16BitSample(smp: Double): SmallInt;
begin
  Result := EnsureRange(round(smp), Low(SmallInt), High(SmallInt));
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
      writeln(ExtractFileName(ParamStr(0)) + ' <source file> <dest file> [quality 0.0-1.0] [iter count 1-inf] [min chunk size 1-inf]');
      WriteLn;
      Exit;
    end;

    enc := TEncoder.Create;

    enc.quality := EnsureRange(StrToFloatDef(ParamStr(3), 0.5), 0.001, 1.0);
    enc.restartCount := StrToIntDef(ParamStr(4), 10);
    enc.minChunkSize := StrToIntDef(ParamStr(5), 32);

    try

      enc.Load(ParamStr(1));
      enc.MakeBands;
      enc.MakeDstData;
      enc.Save(ParamStr(2));

    finally
      enc.Free;
    end;

{$if true}
    ShellExecute(0, 'open', PAnsiChar(ParamStr(2)), nil, nil, 0);
{$else}
    DoExternalEAQUAL(ParamStr(1), ParamStr(2), True, False, 2048);
{$endif}

    ReadLn;
  except
    on e: Exception do
    begin
      WriteLn('Exception: ', e.Message, ' (', e.ClassName, ')');
      ReadLn;
    end;
  end;
end.

