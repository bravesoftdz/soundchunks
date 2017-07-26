program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, Types, math, fgl, MTProcs, yakmo;

type

  TEncoder = class;

  { TChunk }

  TChunk = class
  public
    encoder: TEncoder;

    rawData: TSmallIntDynArray;
    reducedChunk: TChunk;
    dct: TDoubleDynArray;
    weightedDct: TDoubleDynArray;

    joinIndex: Integer;
    joinPenalty: Double;
    dctAddCount: Integer;

    constructor Create(enc: TEncoder; idx: Integer);

    procedure ComputeDCT;
    procedure ComputeInvDCT;
    procedure AddToDCT(ch: TChunk);
    procedure FinalizeDCTAdd;

    procedure FindBestJoinIndex(nextChunk: TChunk);
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
    joinSize: Integer;

    srcHeader: array[$00..$2b] of Byte;
    srcData: TSmallIntDynArray;
    srcDataCount: Integer;
    dstData: TSmallIntDynArray;
    chunkList: TChunkList;
    reducedChunks: TChunkList;
    chunkFreqs: TDoubleDynArray;
    chunkWgtAtts: TDoubleDynArray;
    chunkJoinAtts: TDoubleDynArray;

    class function make16BitSample(smp: Double): SmallInt;
    class function ComputeDCT(chunkSz: Integer; const samples: TSmallIntDynArray): TDoubleDynArray;
    class function ComputeInvDCT(chunkSz: Integer; const dct: TDoubleDynArray): TSmallIntDynArray;
    class function CompareDCT(firstCoeff, lastCoeff: Integer; const dctA, dctB: TDoubleDynArray): Double;
    class function CompressDCT(idx: Integer; coeff: Double): Double;

    constructor Create;
    destructor Destroy; override;

    procedure kMeansReduce;
    procedure load(fn: String);
    procedure makeChunks;
    procedure makeDstData;
    procedure save(fn: String);
  end;


{ TChunk }

constructor TChunk.Create(enc: TEncoder; idx: Integer);
begin
  encoder := enc;
  SetLength(rawData, encoder.chunkSize);
  Move(encoder.srcData[idx * (encoder.chunkSize - encoder.joinSize)], rawData[0], encoder.chunkSize * SizeOf(SmallInt));
  reducedChunk := Self;
  SetLength(dct, encoder.chunkSize);
  SetLength(weightedDct, encoder.chunkSize);
end;

procedure TChunk.ComputeDCT;
var
  k: Integer;
  dta: TSmallIntDynArray;
begin
  dta := Copy(rawData);

  //for k := 0 to encoder.chunkSize - 1 do
  //  dta[k] := TEncoder.make16BitSample(dta[k] * encoder.chunkJoinAtts[k]);

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

procedure TChunk.FindBestJoinIndex(nextChunk: TChunk);
var
  i, j, k, l: Integer;
  smps: TSmallIntDynArray;
  dctOri, dctIter: TDoubleDynArray;
  cmp, curAtt, nextAtt: Double;
begin
  SetLength(smps, encoder.joinSize);
  for i := 0 to encoder.joinSize - 1 do
    smps[i] := rawData[encoder.chunkSize - encoder.joinSize + i];
  dctOri := TEncoder.ComputeDCT(encoder.joinSize, smps);

  // find the best splice point

  joinIndex := 0;
  joinPenalty := MaxDouble;

  i := encoder.joinSize div 2;
  for i := 0 to encoder.joinSize - 1 do
  begin
    l := i - encoder.joinSize div 2;

    for j := 0 to encoder.joinSize - 1 do
    begin
      k := encoder.chunkSize - encoder.joinSize + j;
      curAtt := encoder.chunkJoinAtts[EnsureRange(k + l, encoder.chunkSize - encoder.joinSize, encoder.chunkSize - 1)];
      nextAtt := encoder.chunkJoinAtts[EnsureRange(j + l, 0, encoder.joinSize - 1)];

      Assert(SameValue(curAtt + nextAtt, 1.0));

      smps[j] := TEncoder.make16BitSample(curAtt * reducedChunk.rawData[k] + nextAtt * nextChunk.reducedChunk.rawData[j]);
    end;

    dctIter := TEncoder.ComputeDCT(encoder.joinSize, smps);

    cmp := TEncoder.CompareDCT(0, encoder.joinSize - 1, dctOri, dctIter);
    if cmp < joinPenalty then
    begin
      joinPenalty := cmp;
      joinIndex := i;
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
    //chunkList[AIndex].dct[0] := 0; // remove DC
    //chunk.ComputeInvDCT;
  end;
var
  i, att: Integer;
  chunk: TChunk;
  chunkCount: Integer;
  fsq: Double;
begin
  chunkCount := (srcDataCount - 1) div (chunkSize - joinSize) + 1;
  desiredChunkCount := round(chunkCount * quality);

  WriteLn('makeChunks ', chunkSize, ' * ', chunkCount);

  SetLength(chunkJoinAtts, chunkSize);
  SetLength(chunkWgtAtts, chunkSize);
  SetLength(chunkFreqs, chunkSize);
  for i := 0 to chunkSize - 1 do
  begin
    chunkFreqs[i] := i * sampleRate / (chunkSize - 1) / 4.0;
    fsq := sqr(max(chunkFreqs[i], 20.0));
    chunkWgtAtts[i] := sqr(12194.0) * fsq / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)));

    if i <= joinSize - 1 then
      att := joinSize - i - 1
    else if i >= chunkSize - joinSize - 1 then
      att := i - (chunkSize - joinSize - 1)
    else
      att := 0;

    if joinSize = 0 then
    begin
      chunkJoinAtts[i] := 1.0;
    end
    else
    begin
      chunkJoinAtts[i] := (joinSize - att) / joinSize;
    end;

    //writeln(FloatToStr(att), #9, FloatToStr(chunkJoinAtts[i]));
  end;

  chunkList.Capacity := chunkCount;
  for i := 0 to chunkCount - 1 do
  begin
    chunk := TChunk.Create(Self, i);
    chunkList.Add(chunk);
  end;

  ProcThreadPool.DoParallelLocalProc(@DoDCT, 0, chunkList.Count - 1, nil);
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

  FN := ExtractFilePath(ParamStr(0)) + 'dataset-'+IntToStr(GetCurrentThreadId)+'.txt';
  Dataset := TStringList.Create;
  Dataset.LineBreak := #10;

  try
    for i := 0 to chunkList.Count - 1 do
    begin
      Line := IntToStr(i) + ' ';
      for j := 0 to chunkSize - 1 do
      begin
        v1 := chunkList[i].dct[j];
        v1 := CompressDCT(j, v1);
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
    chunkList[AIndex].FindBestJoinIndex(chunkList[AIndex + 1]);
  end;

var
  i, j, k, l, pos: Integer;
  smp, curAtt, nextAtt: Double;
  penalty: Double;
begin
  WriteLn('makeDstData');

  ProcThreadPool.DoParallelLocalProc(@DoFind, 0, chunkList.Count - 2, nil);

  penalty := 0.0;
  SetLength(dstData, Length(srcData));
  FillWord(dstData[0], Length(srcData), 0);
  for i := 0 to chunkList.Count - 2 do
  begin
    l := chunkList[i].joinIndex - joinSize div 2;
    penalty += chunkList[i].joinPenalty;

    for k := joinSize to chunkSize - 1 do
    begin
      pos := i * (chunkSize - joinSize) + k;
      j := k - chunkSize + joinSize;

      if j >= 0 then
      begin
        curAtt := chunkJoinAtts[EnsureRange(k + l, chunkSize - joinSize, chunkSize - 1)];
        nextAtt := chunkJoinAtts[EnsureRange(j + l, 0, joinSize - 1)];

        Assert(SameValue(curAtt + nextAtt, 1.0));

        smp := curAtt * chunkList[i].reducedChunk.rawData[k] + nextAtt * chunkList[i + 1].reducedChunk.rawData[j];
      end
      else
      begin
        smp := chunkList[i].reducedChunk.rawData[k];
      end;

      dstData[pos] := TEncoder.make16BitSample(smp);
    end;
  end;

  WriteLn(FloatToStr(penalty / (chunkList.Count - 1)), ' average join penalty');
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

class function TEncoder.CompareDCT(firstCoeff, lastCoeff: Integer; const dctA, dctB: TDoubleDynArray): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := firstCoeff to lastCoeff do
    Result += sqr(CompressDCT(i, dctA[i]) - CompressDCT(i, dctB[i]));
  Result := sqrt(Result);
end;

class function TEncoder.CompressDCT(idx: Integer; coeff: Double): Double;
begin
  Result := Sign(coeff) * power(Abs(coeff), 0.707);
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
      writeln(ExtractFileName(ParamStr(0)) + ' <source file> <dest file> [chunk size] [quality 0.0-1.0] [kmeans restarts 1-inf] [join ratio 0.0-0.5]');
      WriteLn;
      Exit;
    end;

    enc := TEncoder.Create;

    enc.quality := EnsureRange(StrToFloatDef(ParamStr(3), 0.5), 0.001, 1.0);
    enc.chunkSize := StrToIntDef(ParamStr(4), 256);
    enc.restartCount := StrToIntDef(ParamStr(5), 4);
    enc.joinSize := EnsureRange(round(StrToFloatDef(ParamStr(6), 0.125) * enc.chunkSize), 0, enc.chunkSize div 2);

    try

      enc.load(ParamStr(1));

      enc.makeChunks;

      enc.kMeansReduce;

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

