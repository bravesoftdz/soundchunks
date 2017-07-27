program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, Types, fgl, MTProcs, math, yakmo, fft, ap, conv;

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
    srcModulo: Integer;
    srcData: TSmallIntDynArray;
    reducedChunk: TChunk;

    fft: TComplex1DArray;
    weightedFft: TComplex1DArray;

    joinPhase: Integer;
    joinPenalty: Double;
    fftAddCount: Integer;

    constructor Create(enc: TEncoder; idx: Integer);

    procedure ComputeFFT;
    procedure AddToFFT(ch: TChunk);
    procedure FinalizeFFTAdd;

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
    highPassData: TSmallIntDynArray;
    chunkList: TChunkList;
    reducedChunks: TChunkList;
    chunkFreqs: TDoubleDynArray;
    chunkWgtAtts: TDoubleDynArray;

    class function make16BitSample(smp: Double): SmallInt;
    class function CompareFFT(firstCoeff, lastCoeff: Integer; compress: Boolean; const fftA, fftB: TDoubleDynArray): Double;
    class function CompressFFT(coeff: Double): Double;
    class function CheckJoinPenalty(x, y, z, a, b, c: Double; TestRange: Boolean): Boolean; inline;
    class function ComputeFFT(chunkSz: Integer; const samples: TSmallIntDynArray): TComplex1DArray;
    class function ComputeInvFFT(chunkSz: Integer; const fft: TComplex1DArray): TSmallIntDynArray;

    constructor Create;
    destructor Destroy; override;

    procedure kMeansReduce;
    procedure load(fn: String);
    procedure makeChunks;
    procedure makeDstData;
    procedure save(fn: String);

    function DoHPFilter(chunkSz: Integer; const samples: TSmallIntDynArray): TSmallIntDynArray;

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
  srcModulo := encoder.chunkSize;

  SetLength(srcData, encoder.chunkSize);
  Move(encoder.highPassData[idx * encoder.chunkSize], srcData[0], encoder.chunkSize * SizeOf(SmallInt));

  reducedChunk := Self;
  SetLength(fft, encoder.chunkSize);
  SetLength(weightedFft, encoder.chunkSize);
end;

procedure TChunk.ComputeFFT;
var
  k: Integer;
  dta: TSmallIntDynArray;
begin
  SetLength(dta, encoder.chunkSize);

  for k := 0 to encoder.chunkSize - 1 do
    dta[k] := srcData[k mod srcModulo];

  fft := encoder.ComputeFFT(encoder.chunkSize, dta);

  for k := 0 to encoder.chunkSize - 1 do
  begin
    weightedFft[k] := fft[k];
    weightedFft[k].X *= encoder.chunkWgtAtts[k];
  end;
end;

procedure TChunk.AddToFFT(ch: TChunk);
var
  k: Integer;
  lfft: TComplex1DArray;
begin
  lfft := encoder.ComputeFFT(encoder.chunkSize, ch.srcData);
  for k := 0 to encoder.chunkSize - 1 do
  begin
    fft[k].X += lfft[k].X;
    fft[k].Y += lfft[k].Y;
  end;
  Inc(fftAddCount);
end;

procedure TChunk.FinalizeFFTAdd;
var
  k: Integer;
begin
  if fftAddCount = 0 then Exit;

  for k := 0 to encoder.chunkSize - 1 do
  begin
    fft[k].X /= fftAddCount;
    fft[k].Y /= fftAddCount;
  end;

  fftAddCount := 0;

  srcData := encoder.ComputeInvFFT(encoder.chunkSize, fft);

  FindBestModulo;
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

    x := srcData[0];
    y := srcData[PhaseSearch];
    z := srcData[PhaseSearch * 2];

    a := srcData[i];
    b := srcData[i + PhaseSearch];
    c := srcData[i + PhaseSearch * 2];

    if not TEncoder.CheckJoinPenalty(x, y, z, a, b, c, True) then
      Continue;

    for j := 0 to encoder.chunkSize - 1 do
      smps[cnt, j] := srcData[j mod i];

    mods[cnt] := i;

    Inc(cnt);
  end;

  SetLength(smps, cnt);
  SetLength(mods, cnt);

  res := encoder.ComputeEAQUALMulti(encoder.chunkSize, False, srcData, smps);

  srcModulo := encoder.chunkSize;
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

  pmd := prevChunk.reducedChunk.srcModulo;
  md := reducedChunk.srcModulo;

  SetLength(smpOri, encoder.chunkSize * 2);
  SetLength(smpItr, encoder.chunkSize, encoder.chunkSize * 2);
  SetLength(phs, encoder.chunkSize);

  move(prevChunk.srcData[0], smpOri[0], encoder.chunkSize * SizeOf(SmallInt));
  move(srcData[0], smpOri[encoder.chunkSize], encoder.chunkSize * SizeOf(SmallInt));

  for i := 0 to encoder.chunkSize - 1 do
    for j := 0 to encoder.chunkSize - 1 do
      smpItr[i, j] := prevChunk.reducedChunk.srcData[j mod pmd];

  checkJP := True;
  repeat
    cnt := 0;

    for i := 0 to md - 1 do
    begin
      x := prevChunk.reducedChunk.srcData[(encoder.chunkSize) mod pmd];
      y := prevChunk.reducedChunk.srcData[(encoder.chunkSize + PhaseSearch) mod pmd];
      z := prevChunk.reducedChunk.srcData[(encoder.chunkSize + PhaseSearch * 2) mod pmd];

      a := reducedChunk.srcData[i];
      b := reducedChunk.srcData[(i + PhaseSearch) mod md];
      c := reducedChunk.srcData[(i + PhaseSearch * 2) mod md];

      if (i <> 0) and not TEncoder.CheckJoinPenalty(x, y, z, a, b, c, checkJP) then
        Continue;

      for j := 0 to encoder.chunkSize - 1 do
        smpItr[cnt, j + encoder.chunkSize] := reducedChunk.srcData[(j + i) mod md];

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

  procedure DoFFT(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  begin
    chunkList[AIndex].ComputeFFT;
    chunkList[AIndex].FindBestModulo;
  end;
var
  i: Integer;
  chunk: TChunk;
  chunkCount: Integer;
  fsq: Double;
begin
  chunkCount := (srcDataCount - 1) div chunkSize + 1;
  desiredChunkCount := round(chunkCount * quality);

  highPassData := DoHPFilter(chunkSize, srcData);

  WriteLn('makeChunks ', chunkSize, ' * ', chunkCount);

  SetLength(chunkWgtAtts, chunkSize);
  SetLength(chunkFreqs, chunkSize);
  for i := 0 to chunkSize - 1 do
  begin
    chunkFreqs[i] := i * sampleRate / (chunkSize - 1) / 2.0;
    fsq := system.sqr(max(chunkFreqs[i], 20.0));
    chunkWgtAtts[i] := system.sqr(12194.0) * fsq / ((fsq + system.sqr(20.6)) * (fsq + system.sqr(12194.0)));
  end;

  chunkList.Capacity := chunkCount;
  for i := 0 to chunkCount - 1 do
  begin
    chunk := TChunk.Create(Self, i);
    chunkList.Add(chunk);
  end;

  ProcThreadPool.DoParallelLocalProc(@DoFFT, 0, chunkList.Count - 1, nil);
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
          reducedChunk.AddToFFT(chunkList[i]);
          chunkList[i].reducedChunk := reducedChunk;
        end;

    reducedChunk.FinalizeFFTAdd;
  end;

var
  FN, Line: String;
  v1, v2: Double;
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
        v1 := CompressFFT(chunkList[i].fft[j].X);
        v2 := CompressFFT(chunkList[i].fft[j].Y);
        Line := Line + Format('%d:%.12g %d:%.12g ', [j * 2, v1, j * 2 + 1, v2]);
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

    writeln(chunk.reducedChunk.srcModulo, #9, chunk.joinPhase, #9, phase);

    for j := 0 to chunkSize - 1 do
      dstData[i * chunkSize + j] := chunk.reducedChunk.srcData[(phase + chunk.joinPhase + j) mod chunk.reducedChunk.srcModulo];

    //phase := (phase + chunk.joinPhase) mod chunk.reducedChunk.srcModulo;
  end;

  WriteLn('avg join penalty: ', FloatToStr(penalty / chunkList.Count));
end;

class function TEncoder.make16BitSample(smp: Double): SmallInt;
begin
  Result := EnsureRange(round(smp), Low(SmallInt), High(SmallInt));
end;

class function TEncoder.ComputeFFT(chunkSz: Integer; const samples: TSmallIntDynArray): TComplex1DArray;
var
  k: Integer;
  rl: TReal1DArray;
begin
  SetLength(rl, chunkSz);

  for k := 0 to chunkSz - 1 do
    rl[k] := samples[k];

  FFTR1D(rl, chunkSz, Result);
end;

class function TEncoder.ComputeInvFFT(chunkSz: Integer; const fft: TComplex1DArray): TSmallIntDynArray;
var
  k: Integer;
  ifft: TDoubleDynArray;
begin
  FFTR1DInv(fft, chunkSz, ifft);

  SetLength(Result, chunkSz);
  for k := 0 to chunkSz - 1 do
    Result[k] := make16BitSample(ifft[k]);
end;

function TEncoder.DoHPFilter(chunkSz: Integer; const samples: TSmallIntDynArray): TSmallIntDynArray;
var
  fc, b, sinc, win, sum: Double;
  i, N: Integer;
  h: TReal1DArray;
  ins, res: TReal1DArray;
begin
  fc := 1.0 / chunkSz;
  b := fc * 0.001;

  fc += b;
  N := ceil(4 / b);
  if (N mod 2) = 0 then N += 1;
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

  for i := 0 to N - 1 do
    h[i] := -h[i] / sum;

  h[(N - 1) div 2] += 1;

  writeln('DoHPFilter ', FloatToStr(sampleRate * fc), ' ', N);

  SetLength(ins, Length(samples));
  for i := 0 to High(samples) do
    ins[i] := samples[i];

  ConvR1D(ins, Length(samples), h, N, res);

  SetLength(Result, Length(samples));
  for i := 0 to High(samples) do
    Result[i] := make16BitSample(res[i + (N - 1) div 2]);
end;

class function TEncoder.CompareFFT(firstCoeff, lastCoeff: Integer; compress: Boolean; const fftA, fftB: TDoubleDynArray): Double;
var
  i: Integer;
begin
  Result := 0.0;
  if compress then
  begin
    for i := firstCoeff to lastCoeff do
      Result += system.sqr(CompressFFT(fftA[i]) - CompressFFT(fftB[i]));
  end
  else
  begin
    for i := firstCoeff to lastCoeff do
      Result += system.sqr(fftA[i] - fftB[i]);
  end;
  Result := system.sqrt(Result);
end;

class function TEncoder.CompressFFT(coeff: Double): Double;
begin
  Result := Sign(coeff) * system.sqrt(system.Abs(coeff));
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

