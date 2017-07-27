program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, Types, math, fgl, MTProcs, yakmo;

const
  PhaseSearch = 1;

type
  TSmallIntDynArray2 = array of TSmallIntDynArray;

  TEncoder = class;

  { TChunk }

  TChunk = class
  public
    encoder: TEncoder;

    index: Integer;
    rawModulo: Integer;
    rawData: TSmallIntDynArray;
    reducedChunk: TChunk;
    dct: TDoubleDynArray;
    weightedDct: TDoubleDynArray;

    joinPhase: Integer;
    joinPenalty: Double;
    dctAddCount: Integer;

    constructor Create(enc: TEncoder; idx: Integer);

    procedure ComputeDCT;
    procedure ComputeInvDCT;
    procedure AddToDCT(ch: TChunk);
    procedure FinalizeDCTAdd;

    procedure FindBestModulo;
    procedure FindBestJoinPhase(prevChunk: TChunk);
  end;

  TChunkList = specialize TFPGObjectList<TChunk>;

  { TEncoder }

  TEncoder = class
  public
    chunkSize: Integer;
    quality: Double;
    desiredChunkCount: Integer;
    restartCount: Integer;
    sampleRate: Integer;

    srcHeader: array[$00..$2b] of Byte;
    srcData: TSmallIntDynArray;
    srcDataCount: Integer;
    dstData: TSmallIntDynArray;
    chunkList: TChunkList;
    reducedChunks: TChunkList;
    chunkFreqs: TDoubleDynArray;
    chunkWgtAtts: TDoubleDynArray;

    class function make16BitSample(smp: Double): SmallInt;
    class function ComputeDCT(chunkSz: Integer; const samples: TSmallIntDynArray): TDoubleDynArray;
    class function ComputeInvDCT(chunkSz: Integer; const dct: TDoubleDynArray): TSmallIntDynArray;
    class function CompareDCT(firstCoeff, lastCoeff: Integer; compress: Boolean; const dctA, dctB: TDoubleDynArray): Double;
    class function CompressDCT(coeff: Double): Double;
    class function CheckJoinPenalty(x, y, z, a, b, c: Double; TestRange: Boolean): Boolean; inline;

    constructor Create;
    destructor Destroy; override;

    procedure kMeansReduce;
    procedure load(fn: String);
    procedure makeChunks;
    procedure makeDstData;
    procedure save(fn: String);

    function ComputeEAQUAL(chunkSz: Integer; UseDIX: Boolean; const smpRef, smpTst: TSmallIntDynArray): Double;
    function ComputeEAQUALMulti(chunkSz: Integer; UseDIX: Boolean; const smpRef: TSmallIntDynArray;
      smpTst: TSmallIntDynArray2): TDoubleDynArray;
  end;

function Div0(x, y: Double): Double; inline;
begin
  Result := 0.0;
  if not IsZero(y) then
    Result := x / y;
end;

{ TChunk }

constructor TChunk.Create(enc: TEncoder; idx: Integer);
begin
  index := idx;
  encoder := enc;
  rawModulo := encoder.chunkSize;
  SetLength(rawData, encoder.chunkSize);
  Move(encoder.srcData[idx * encoder.chunkSize], rawData[0], encoder.chunkSize * SizeOf(SmallInt));
  reducedChunk := Self;
  SetLength(dct, encoder.chunkSize);
  SetLength(weightedDct, encoder.chunkSize);
end;

procedure TChunk.ComputeDCT;
var
  k: Integer;
  dta: TSmallIntDynArray;
begin
  SetLength(dta, encoder.chunkSize);

  for k := 0 to encoder.chunkSize - 1 do
    dta[k] := TEncoder.make16BitSample(rawData[k mod rawModulo]);

  dct := TEncoder.ComputeDCT(encoder.chunkSize, dta);

  for k := 0 to encoder.chunkSize - 1 do
    weightedDct[k] := dct[k] * encoder.chunkWgtAtts[k];
end;

procedure TChunk.ComputeInvDCT;
begin
  rawData := TEncoder.ComputeInvDCT(encoder.chunkSize, dct);
end;

procedure TChunk.AddToDCT(ch: TChunk);
var
  k: Integer;
  ldct: TDoubleDynArray;
begin
  ldct := TEncoder.ComputeDCT(encoder.chunkSize, ch.rawData);
  for k := 0 to encoder.chunkSize - 1 do
    dct[k] += ldct[k];
  Inc(dctAddCount);
end;

procedure TChunk.FinalizeDCTAdd;
var
  k: Integer;
begin
  if dctAddCount = 0 then Exit;

  for k := 0 to encoder.chunkSize - 1 do
    dct[k] /= dctAddCount;

  dctAddCount := 0;

  ComputeInvDCT;
end;

procedure TChunk.FindBestModulo;
var
  i, j, a, b, c, x, y, z, cnt: Integer;
  cmp, best: Double;
  smps: TSmallIntDynArray2;
  mods: TIntegerDynArray;
  res: TDoubleDynArray;
  checkJP: Boolean;
begin
  // find the best looping point in a chunk

  if index and $0f = 0 then Write('.');

  SetLength(smps, encoder.chunkSize, encoder.chunkSize);
  SetLength(mods, encoder.chunkSize);
  cnt := 0;

  for i := encoder.chunkSize - PhaseSearch * 2 - 1 downto encoder.chunkSize div 8 do
  begin

    x := rawData[0];
    y := rawData[PhaseSearch];
    z := rawData[PhaseSearch * 2];

    a := rawData[i];
    b := rawData[i + PhaseSearch];
    c := rawData[i + PhaseSearch * 2];

    if not TEncoder.CheckJoinPenalty(x, y, z, a, b, c, True) then
      Continue;

    for j := 0 to encoder.chunkSize - 1 do
      smps[cnt, j] := rawData[j mod i];

    mods[cnt] := i;

    Inc(cnt);
  end;

  SetLength(smps, cnt);
  SetLength(mods, cnt);

  res := encoder.ComputeEAQUALMulti(encoder.chunkSize, False, rawData, smps);

  rawModulo := encoder.chunkSize;
  best := MaxDouble;

  for i := 0 to cnt - 1 do
  begin
    cmp := res[i];

    if cmp < best then
    begin
      best := cmp;
      rawModulo := mods[i];
    end;
  end;
end;

procedure TChunk.FindBestJoinPhase(prevChunk: TChunk);
var
  i, j, a, b, c, x, y, z, pmd, md, cnt: Integer;
  cmp: Double;
  checkJP: Boolean;
  smpOri: TSmallIntDynArray;
  smpItr: TSmallIntDynArray2;
  phs: TIntegerDynArray;
  res: TDoubleDynArray;
begin
  // find the best phase shift to join the chunks

  if index and $0f = 0 then Write('.');

  pmd := prevChunk.reducedChunk.rawModulo;
  md := reducedChunk.rawModulo;

  SetLength(smpOri, encoder.chunkSize * 2);
  SetLength(smpItr, encoder.chunkSize, encoder.chunkSize * 2);
  SetLength(phs, encoder.chunkSize);

  move(prevChunk.rawData[0], smpOri[0], encoder.chunkSize * SizeOf(SmallInt));
  move(rawData[0], smpOri[encoder.chunkSize], encoder.chunkSize * SizeOf(SmallInt));

  for i := 0 to encoder.chunkSize - 1 do
    for j := 0 to encoder.chunkSize - 1 do
      smpItr[i, j] := prevChunk.reducedChunk.rawData[j mod pmd];

  checkJP := True;
  repeat
    cnt := 0;

    for i := 0 to md - 1 do
    begin
      x := prevChunk.reducedChunk.rawData[(encoder.chunkSize) mod pmd];
      y := prevChunk.reducedChunk.rawData[(encoder.chunkSize + PhaseSearch) mod pmd];
      z := prevChunk.reducedChunk.rawData[(encoder.chunkSize + PhaseSearch * 2) mod pmd];

      a := reducedChunk.rawData[i];
      b := reducedChunk.rawData[(i + PhaseSearch) mod md];
      c := reducedChunk.rawData[(i + PhaseSearch * 2) mod md];

      if (i <> 0) and not TEncoder.CheckJoinPenalty(x, y, z, a, b, c, checkJP) then
        Continue;

      for j := 0 to encoder.chunkSize - 1 do
        smpItr[cnt, j + encoder.chunkSize] := reducedChunk.rawData[(j + i) mod md];

      phs[cnt] := i;
      Inc(cnt);
    end;

    checkJP := False;
  until cnt <> 0;

  SetLength(smpItr, cnt);
  SetLength(phs, cnt);

  res := encoder.ComputeEAQUALMulti(encoder.chunkSize * 2, False, smpOri, smpItr);

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

procedure TEncoder.load(fn: String);
var
  fs: TFileStream;
begin
  WriteLn('load ', fn);
  fs := TFileStream.Create(fn, fmOpenRead or fmShareDenyNone);
  try
    fs.ReadBuffer(srcHeader[0], SizeOf(srcHeader));
    srcDataCount := (fs.Size - fs.Position) div 2;
    SetLength(srcData, srcDataCount + chunkSize);
    FillWord(srcData[0], srcDataCount + chunkSize, 0);

    fs.ReadBuffer(srcData[0], srcDataCount * 2);
  finally
    fs.Free;
  end;

  sampleRate := PInteger(@srcHeader[$18])^;
  writeln(sampleRate, ' Hz');
end;

procedure TEncoder.save(fn: String);
var
  fs: TFileStream;
begin
  WriteLn('save ', fn);

  fs := TFileStream.Create(fn, fmCreate or fmShareDenyWrite);
  try
    fs.WriteBuffer(srcHeader[0], SizeOf(srcHeader));
    fs.WriteBuffer(dstData[0], srcDataCount * 2);
  finally
    fs.Free;
  end;
end;

constructor TEncoder.Create;
begin
  chunkList := TChunkList.Create;
  reducedChunks := TChunkList.Create;
end;

destructor TEncoder.Destroy;
begin
  chunkList.Free;
  reducedChunks.Free;

  inherited Destroy;
end;

procedure TEncoder.makeChunks;

  procedure DoDCT(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    chunk: TChunk;
  begin
    chunk:= chunkList[AIndex];

    chunk.ComputeDCT;
    //chunk.dct[0] := 0;
    //chunk.weightedDct[0] := 0;
    //chunk.ComputeInvDCT;
    chunk.FindBestModulo;
  end;
var
  i: Integer;
  chunk: TChunk;
  chunkCount: Integer;
  fsq: Double;
begin
  chunkCount := (srcDataCount - 1) div chunkSize + 1;
  desiredChunkCount := round(chunkCount * quality);

  WriteLn('makeChunks ', chunkSize, ' * ', chunkCount);

  SetLength(chunkWgtAtts, chunkSize);
  SetLength(chunkFreqs, chunkSize);
  for i := 0 to chunkSize - 1 do
  begin
    chunkFreqs[i] := i * sampleRate / (chunkSize - 1) / 4.0;
    fsq := sqr(max(chunkFreqs[i], 20.0));
    chunkWgtAtts[i] := sqr(12194.0) * fsq / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)));
  end;

  chunkList.Capacity := chunkCount;
  for i := 0 to chunkCount - 1 do
  begin
    chunk := TChunk.Create(Self, i);
    chunkList.Add(chunk);
  end;

  ProcThreadPool.DoParallelLocalProc(@DoDCT, 0, chunkList.Count - 1, nil);
  WriteLn;
end;

procedure TEncoder.kMeansReduce;
var
  XYC: TIntegerDynArray;

  procedure DoXYC(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
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
  WriteLn('kMeansReduce ', desiredChunkCount);

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
        v1 := CompressDCT(v1);
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
  DoExternalKMeans(FN, desiredChunkCount, RestartCount, False, XYC);

  for i := 0 to desiredChunkCount - 1 do
    reducedChunks.Add(TChunk.Create(Self, 0));

  ProcThreadPool.DoParallelLocalProc(@DoXYC, 0, desiredChunkCount - 1, nil);
end;

procedure TEncoder.makeDstData;

  procedure DoFind(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  begin
    chunkList[AIndex].FindBestJoinPhase(chunkList[AIndex - 1]);
  end;

var
  i, j: Integer;
  penalty: Double;
  chunk: TChunk;
  phase: Integer;
begin
  WriteLn('makeDstData');

  //for i := 1 to chunkList.Count - 1 do
  //  chunkList[i].FindBestJoinPhase(chunkList[i - 1]);
  //
  ProcThreadPool.DoParallelLocalProc(@DoFind, 1, chunkList.Count - 1, nil);
  WriteLn;

  phase := 0;
  penalty := 0.0;
  SetLength(dstData, Length(srcData));
  FillWord(dstData[0], Length(srcData), 0);
  for i := 0 to chunkList.Count - 1 do
  begin
    chunk := chunkList[i];

    if not IsInfinite(chunk.joinPenalty) then
      penalty += chunk.joinPenalty;

    writeln(chunk.reducedChunk.rawModulo, #9, chunk.joinPhase, #9, phase);

    for j := 0 to chunkSize - 1 do
      dstData[i * chunkSize + j] := chunk.reducedChunk.rawData[(phase + chunk.joinPhase + j) mod chunk.reducedChunk.rawModulo];

//    phase += chunk.joinPhase;
  end;

  WriteLn('avg join penalty: ', FloatToStr(penalty / chunkList.Count));
end;

class function TEncoder.make16BitSample(smp: Double): SmallInt;
begin
  Result := EnsureRange(round(smp), Low(SmallInt), High(SmallInt));
end;

class function TEncoder.ComputeDCT(chunkSz: Integer; const samples: TSmallIntDynArray): TDoubleDynArray;
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

class function TEncoder.ComputeInvDCT(chunkSz: Integer; const dct: TDoubleDynArray): TSmallIntDynArray;
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
    Result[n] := TEncoder.make16BitSample(sum * sqrt(2.0 / chunkSz));
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

function TEncoder.ComputeEAQUAL(chunkSz: Integer; UseDIX: Boolean; const smpRef, smpTst: TSmallIntDynArray): Double;
var
  FNRef, FNTst: String;
  ms: TMemoryStream;
begin

  FNRef := GetTempFileName('', 'ref-'+IntToStr(GetCurrentThreadId)+'.wav');
  FNTst := GetTempFileName('', 'tst-'+IntToStr(GetCurrentThreadId)+'.wav');

  ms := TMemoryStream.Create;
  try
    ms.Write(srcHeader[0], $28);
    ms.WriteDWord(chunkSz * SizeOf(SmallInt));
    ms.Write(smpRef[0], chunkSz * SizeOf(SmallInt));

    ms.SaveToFile(FNRef);
    ms.Clear;

    ms.Write(srcHeader[0], $28);
    ms.WriteDWord(chunkSz * SizeOf(SmallInt));
    ms.Write(smpTst[0], chunkSz * SizeOf(SmallInt));

    ms.SaveToFile(FNTst);
  finally
    ms.Free;
  end;

  Result := DoExternalEAQUAL(FNRef, FNTst, UseDIX, chunkSz * 2);

  DeleteFile(FNRef);
  DeleteFile(FNTst);
end;

function TEncoder.ComputeEAQUALMulti(chunkSz: Integer; UseDIX: Boolean; const smpRef: TSmallIntDynArray; smpTst: TSmallIntDynArray2
  ): TDoubleDynArray;
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
      ms.Write(smpRef[0], chunkSz * SizeOf(SmallInt));
      ms.Write(zeroes[0], chunkSz * SizeOf(SmallInt));
    end;

    ms.SaveToFile(FNRef);
    ms.Clear;

    ms.Write(srcHeader[0], $28);
    ms.WriteDWord(chunkSz * 2 * SizeOf(SmallInt) * Length(smpTst));
    for i := 0 to High(smpTst) do
    begin
      ms.Write(smpTst[i, 0], chunkSz * SizeOf(SmallInt));
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
      writeln(ExtractFileName(ParamStr(0)) + ' <source file> <dest file> [quality 0.0-1.0] [chunk size] [iter count 1-inf]');
      WriteLn;
      Exit;
    end;

    enc := TEncoder.Create;

    enc.quality := EnsureRange(StrToFloatDef(ParamStr(3), 0.5), 0.001, 1.0);
    enc.chunkSize := EnsureRange(StrToIntDef(ParamStr(4), 8), 16, 4096);
    enc.restartCount := StrToIntDef(ParamStr(5), 4);

    try

      enc.load(ParamStr(1));

      enc.makeChunks;

      //enc.kMeansReduce;

      enc.makeDstData;

      enc.save(ParamStr(2));

    finally
      enc.Free;
    end;

    ShellExecute(0, 'open', PAnsiChar(ParamStr(2)), nil, nil, 0);
    Sleep(5000);
  except
    on e: Exception do
    begin
      WriteLn('Exception: ', e.Message, ' (', e.ClassName, ')');
      ReadLn;
    end;
  end;
end.

