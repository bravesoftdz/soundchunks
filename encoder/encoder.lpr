program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, yakmo, ap, conv;

const
  PhaseSearch = 1;
  BandCount = 4;
  LowCut = 30.0;
  HighCut = 18000.0;

type
  TSmallIntDynArray2 = array of TSmallIntDynArray;
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
    srcModulo: Integer;

    srcData: TDoubleDynArray;
    dct: TDoubleDynArray;

    joinPhase: Integer;
    joinPenalty: Double;
    dctAddCount: Integer;

    constructor Create(enc: TEncoder; bnd: TBand; idx: Integer);

    procedure ComputeDCT;
    procedure InitDCTAdd;
    procedure AddToDCT(ch: TChunk);
    procedure FinalizeDCTAdd;

    procedure FindBestModulo;
    procedure FindBestJoinPhase(prevChunk: TChunk);
  end;

  TChunkList = specialize TFPGObjectList<TChunk>;

  { TBand }

  TBand = class
  public
    encoder: TEncoder;

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

    constructor Create;
    destructor Destroy; override;

    procedure Load(fn: String);
    procedure Save(fn: String);

    procedure FindBestDesiredChunksCounts;
    procedure MakeBands;
    procedure MakeDstData;

    function DoFilter(fc, transFactor: Double; HighPass: Boolean; chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
    function DoBPFilter(fcl, fch, transFactor: Double; chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;

    function ComputeEAQUAL(chunkSz: Integer; UseDIX: Boolean; const smpRef, smpTst: TDoubleDynArray): Double;
    function ComputeEAQUALMulti(chunkSz: Integer; UseDIX: Boolean; const smpRef: TDoubleDynArray;
      smpTst: TDoubleDynArray2): TDoubleDynArray;
  end;

function Div0(x, y: Double): Double; inline;
begin
  Result := 0.0;
  if not IsZero(y) then
    Result := x / y;
end;

constructor TBand.Create(enc: TEncoder; idx: Integer);
var
  i: Integer;
  ratio: Double;
  fcdc: Double;
begin
  encoder := enc;
  index := idx;

  fcdc := LowCut / encoder.sampleRate;
  ratio := log2(HighCut / LowCut) / BandCount;
  fcl := fcdc * power(2.0, index * ratio);
  fch := fcdc * power(2.0, (index + 1) * ratio);

  chunkSize := 1;
  while fcl < 1 / chunkSize do chunkSize *= 2;
  chunkSize := max(encoder.minChunkSize, chunkSize);

  chunkCount := (encoder.srcDataCount - 1) div chunkSize + 1;

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
    //chunk.dct[0] := 0;
    //chunk.srcData := TEncoder.ComputeInvDCT(chunkSize, chunk.dct);

    chunk.FindBestModulo;
  end;

var
  i: Integer;
  chunk: TChunk;
begin
  WriteLn('MakeChunks #', index, ' (', round(fcl * encoder.sampleRate), ' Hz .. ', round(fch * encoder.sampleRate), ' Hz); ', chunkSize, ' * (', chunkCount, ' -> ', desiredChunkCount,')');

  srcData := encoder.DoBPFilter(fcl, fch, 0.001, chunkSize, encoder.srcData);

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
  //exit;

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

procedure TBand.MakeDstData;
var
  i, j: Integer;
  chunk: TChunk;
  phase: Integer;
begin
  WriteLn('MakeDstData #', index);

  for i := 1 to chunkList.Count - 1 do
    chunkList[i].FindBestJoinPhase(chunkList[i - 1]);

  phase := 0;
  SetLength(dstData, Length(srcData));
  FillQWord(dstData[0], Length(srcData), 0);
  for i := 0 to chunkList.Count - 1 do
  begin
    chunk := chunkList[i];

    //writeln(chunk.reducedChunk.srcModulo, #9, chunk.joinPhase, #9, phase);

    for j := 0 to chunkSize - 1 do
      dstData[i * chunkSize + j] := chunk.reducedChunk.srcData[(phase + chunk.joinPhase + j) mod chunk.reducedChunk.srcModulo];

    //phase := (phase + chunk.joinPhase) mod chunk.reducedChunk.srcModulo;
  end;
end;

{ TChunk }

constructor TChunk.Create(enc: TEncoder; bnd: TBand; idx: Integer);
begin
  index := idx;
  encoder := enc;
  band := bnd;
  srcModulo := band.chunkSize;

  SetLength(srcData, band.chunkSize);
  Move(band.srcData[idx * band.chunkSize], srcData[0], band.chunkSize * SizeOf(Double));

  reducedChunk := Self;
  SetLength(dct, band.chunkSize);
end;

procedure TChunk.ComputeDCT;
var
  k: Integer;
  dta: TDoubleDynArray;
begin
  SetLength(dta, band.chunkSize);

  for k := 0 to band.chunkSize - 1 do
    dta[k] := srcData[k mod srcModulo];

  dct := TEncoder.ComputeDCT(band.chunkSize, dta);
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

  FindBestModulo;
end;

procedure TChunk.FindBestModulo;
var
  i, j, cnt: Integer;
  a, b, c, x, y, z: Double;
  cmp, best: Double;
  smps: TDoubleDynArray2;
  mods: TIntegerDynArray;
  res: TDoubleDynArray;
begin
  exit;

  // find the best looping point in a chunk

  if index and $0f = 0 then Write('.');

  SetLength(smps, band.chunkSize, band.chunkSize);
  SetLength(mods, band.chunkSize);
  cnt := 0;

  for i := band.chunkSize - PhaseSearch * 2 - 1 downto band.chunkSize div 8 do
  begin

    x := srcData[0];
    y := srcData[PhaseSearch];
    z := srcData[PhaseSearch * 2];

    a := srcData[i];
    b := srcData[i + PhaseSearch];
    c := srcData[i + PhaseSearch * 2];

    if not TEncoder.CheckJoinPenalty(x, y, z, a, b, c, True) then
      Continue;

    for j := 0 to band.chunkSize - 1 do
      smps[cnt, j] := srcData[j mod i];

    mods[cnt] := i;

    Inc(cnt);
  end;

  SetLength(smps, cnt);
  SetLength(mods, cnt);

  res := encoder.ComputeEAQUALMulti(band.chunkSize, False, srcData, smps);

  srcModulo := band.chunkSize;
  best := MaxDouble;

  for i := 0 to cnt - 1 do
  begin
    cmp := res[i];

    if cmp < best then
    begin
      best := cmp;
      srcModulo := mods[i];
    end;
  end;
end;

procedure TChunk.FindBestJoinPhase(prevChunk: TChunk);
var
  i, j, pmd, md, cnt: Integer;
  a, b, c, x, y, z: Double;
  cmp: Double;
  checkJP: Boolean;
  smpOri: TDoubleDynArray;
  smpItr: TDoubleDynArray2;
  phs: TIntegerDynArray;
  res: TDoubleDynArray;
begin
  exit;

  // find the best phase shift to join the chunks

  if index and $0f = 0 then Write('.');

  pmd := prevChunk.reducedChunk.srcModulo;
  md := reducedChunk.srcModulo;

  SetLength(smpOri, band.chunkSize * 2);
  SetLength(smpItr, band.chunkSize, band.chunkSize * 2);
  SetLength(phs, band.chunkSize);

  move(prevChunk.srcData[0], smpOri[0], band.chunkSize * SizeOf(Double));
  move(srcData[0], smpOri[band.chunkSize], band.chunkSize * SizeOf(Double));

  for i := 0 to band.chunkSize - 1 do
    for j := 0 to band.chunkSize - 1 do
      smpItr[i, j] := prevChunk.reducedChunk.srcData[j mod pmd];

  checkJP := True;
  repeat
    cnt := 0;

    for i := 0 to md - 1 do
    begin
      x := prevChunk.reducedChunk.srcData[(band.chunkSize) mod pmd];
      y := prevChunk.reducedChunk.srcData[(band.chunkSize + PhaseSearch) mod pmd];
      z := prevChunk.reducedChunk.srcData[(band.chunkSize + PhaseSearch * 2) mod pmd];

      a := reducedChunk.srcData[i];
      b := reducedChunk.srcData[(i + PhaseSearch) mod md];
      c := reducedChunk.srcData[(i + PhaseSearch * 2) mod md];

      if (i <> 0) and not TEncoder.CheckJoinPenalty(x, y, z, a, b, c, checkJP) then
        Continue;

      for j := 0 to band.chunkSize - 1 do
        smpItr[cnt, j + band.chunkSize] := reducedChunk.srcData[(j + i) mod md];

      phs[cnt] := i;
      Inc(cnt);
    end;

    checkJP := False;
  until cnt <> 0;

  SetLength(smpItr, cnt);
  SetLength(phs, cnt);

  res := encoder.ComputeEAQUALMulti(band.chunkSize * 2, False, smpOri, smpItr);

  joinPhase := 0;
  joinPenalty := Infinity;

  for i := 0 to cnt - 1 do
  begin
    cmp := res[i];

    if cmp < joinPenalty then
    begin
      joinPenalty := cmp;
      joinPhase := phs[i];
    end;
  end;
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

      bnd.desiredChunkCount := round(sz * dbBatt);
      allSz += bnd.desiredChunkCount * bnd.chunkSize;
    end;
    Inc(sz);
  until allSz >= srcDataCount * quality;

  projectedDataCount := allSz;

  WriteLn('projectedDataCount = ', projectedDataCount);
end;

function TEncoder.DoFilter(fc, transFactor: Double; HighPass: Boolean; chunkSz: Integer; const samples: TDoubleDynArray
  ): TDoubleDynArray;
var
  b, sinc, win, sum: Double;
  i, N: Integer;
  h: TReal1DArray;
begin
  if (fc <= 0.0) and HighPass or (fc >=0.5) and not HighPass then
  begin
    Result := samples;
    Exit;
  end;

  b := fc * transFactor;
  N := ceil(4 / b);
  if (N mod 2) = 0 then N += 1;

  //writeln('DoFilter ', ifthen(HighPass, 'HP', 'LP'), ' ', FloatToStr(sampleRate * fc), ' ', N);

  SetLength(h, N);
  sum := 0;
  for i := 0 to N - 1 do
  begin
    sinc := 2.0 * fc * (i - (N - 1) / 2.0) * pi;
    if IsZero(sinc) then
      sinc := 1.0
    else
      sinc := sin(sinc) / sinc;

    win := 0.42 - 0.5 * cos(2 * pi * i / (N - 1)) + 0.08 * cos(4 * pi * i / (N - 1));

    h[i] := sinc * win;
    sum += h[i];
  end;

  if HighPass then
  begin
    for i := 0 to N - 1 do
      h[i] := -h[i] / sum;

    h[(N - 1) div 2] += 1;
  end
  else
  begin
    for i := 0 to N - 1 do
      h[i] := h[i] / sum;
  end;

  ConvR1D(samples, Length(samples), h, N, Result);

  for i := 0 to High(samples) do
    Result[i] := Result[i + (N - 1) div 2];

  SetLength(Result, Length(samples));
end;

function TEncoder.DoBPFilter(fcl, fch, transFactor: Double; chunkSz: Integer; const samples: TDoubleDynArray
  ): TDoubleDynArray;
begin
  Assert(fch > fcl);
  Result := DoFilter(fcl, transFactor, True, chunkSz, samples);
  Result := DoFilter(fch, transFactor, False, chunkSz, Result);
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
  chunk: TChunk;
  phase: Integer;
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
    enc.minChunkSize := StrToIntDef(ParamStr(5), 64);

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

