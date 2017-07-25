program encoder;

uses Classes, sysutils, Types, math, fgl, MTProcs, yakmo;


type

  { TChunk }

  TChunk = class
  public
    rawData: PSmallInt;
    reducedData: PSmallInt;
    dct: PSingle;
    weightedDct: PSingle;

    constructor Create(idx: Integer);
    destructor Destroy; override;

    procedure ComputeDCT;
  end;

  TChunkList = specialize TFPGObjectList<TChunk>;

var
  sampleRate: Integer;
  srcHeader: array[$00..$2b] of Byte;
  srcData: PSmallInt;
  srcDataCount: Integer;
  dstData: PSmallInt;
  chunkSize: Integer;
  desiredChunkCount: Integer;
  restartCount: Integer;
  chunkList: TChunkList;
  chunkFreqs: PSingle;
  chunkAtts: PSingle;


  { TChunk }

  constructor TChunk.Create(idx: Integer);
  begin
    rawData := @srcData[idx * chunkSize];
    dct := AllocMem(chunkSize * SizeOf(Single));
    weightedDct := AllocMem(chunkSize * SizeOf(Single));
  end;

  destructor TChunk.Destroy;
  begin
    Freemem(dct);
    Freemem(weightedDct);

    inherited Destroy;
  end;

  procedure TChunk.ComputeDCT;
  var
    k, n: Integer;
    sum, s: Single;
  begin
    for k := 0 to chunkSize - 1 do
    begin
      sum := 0;
      s := ifthen(k = 0, sqrt(0.5), 1.0);
      for n := 0 to chunkSize - 1 do
        sum += s * rawData[n] * cos(pi * (n + 0.5) * k / chunkSize);
      dct[k] := sum * sqrt (2.0 / chunkSize);
      weightedDct[k] := dct[k] * chunkAtts[k];
    end;
  end;



procedure load(fn: String);
var
  fs: TFileStream;
  i: Integer;
  prev, cur: SmallInt;
begin
  WriteLn('load ', fn);
  fs := TFileStream.Create(fn, fmOpenRead or fmShareDenyNone);
  try
    fs.ReadBuffer(srcHeader[0], SizeOf(srcHeader));
    srcDataCount := (fs.Size - fs.Position) div 2;
    srcData := AllocMem(srcDataCount * 2);
    fs.ReadBuffer(srcData[0], srcDataCount * 2);
  finally
    fs.Free;
  end;

  prev := 0;
  for i := 0 to srcDataCount - 1 do
  begin
    cur := srcData[i];
    srcData[i] -= prev;
    prev := cur;
  end;

  sampleRate := PInteger(@srcHeader[$18])^;
  writeln(sampleRate, ' Hz');
end;

procedure save(fn: String);
var
  fs: TFileStream;
  i: Integer;
begin
  WriteLn('save ', fn);

  for i := 1 to srcDataCount - 1 do
    dstData[i] += dstData[i - 1];

  fs := TFileStream.Create(fn, fmCreate or fmShareDenyWrite);
  try
    fs.WriteBuffer(srcHeader[0], SizeOf(srcHeader));
    fs.WriteBuffer(dstData[0], srcDataCount * 2);
  finally
    fs.Free;
  end;
end;

procedure makeChunks;

  procedure DoDCT(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  begin
    chunkList[AIndex].ComputeDCT;
  end;
var
  i: Integer;
  chunk: TChunk;
  chunkCount: Integer;
  fsq: Single;
begin
  chunkCount := srcDataCount div chunkSize;

  WriteLn('makeChunks ', chunkSize, ' * ', chunkCount);

  for i := 0 to chunkSize - 1 do
  begin
    chunkFreqs[i] := i * sampleRate / (chunkSize - 1) / 2;
    fsq := sqr(chunkFreqs[i]);
    chunkAtts[i] := sqr(12194.0) * fsq / ((fsq + sqr(20.6)) * (fsq + sqr(12194.0)));
  end;

  chunkList.Capacity := chunkCount;
  for i := 0 to chunkCount - 1 do
  begin
    chunk := TChunk.Create(i);
    chunkList.Add(chunk);
  end;

  ProcThreadPool.DoParallelLocalProc(@DoDCT, 0, chunkList.Count - 1, nil);
end;

procedure kMeansReduce;
var
  FN, Line: String;
  v1: Single;
  Dataset: TStringList;
  XYC: TIntegerDynArray;
  i, j, First: Integer;
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
        v1 := chunkList[i].weightedDct[j];
        Line := Line + Format('%d:%g ', [j, v1]);
      end;
      Dataset.Add(Line);
    end;
    Dataset.SaveToFile(FN);
  finally
    Dataset.Free;
  end;

  SetLength(XYC, chunkList.Count);
  FillChar(XYC[0], chunkList.Count * SizeOF(Integer), $ff);
  DoExternalKMeans(FN, desiredChunkCount, RestartCount, XYC);

  for j := 0 to desiredChunkCount - 1 do
  begin
    First := -1;
    for i := 0 to chunkList.Count - 1 do
      if XYC[i] = j then
      begin
        First := i;
        Break;
      end;

    if First <> -1 then
      for i := 0 to chunkList.Count - 1 do
        if (XYC[i] = j) and (i <> First) then
          chunkList[i].reducedData := chunkList[First].rawData;
  end;
end;

begin
  WriteLn('Usage: (source file must be 16bit mono WAV)');
  writeln(ExtractFileName(ParamStr(0)) + ' <source file> <dest file> <chunk size> <desired chunk count> <kmeans restarts>');

  chunkSize := StrToIntDef(ParamStr(3), 512);
  desiredChunkCount := StrToIntDef(ParamStr(4), 4096);
  restartCount := StrToIntDef(ParamStr(5), 4);

  chunkAtts := AllocMem(chunkSize * SizeOf(Single));
  chunkFreqs := AllocMem(chunkSize * SizeOf(Single));
  chunkList := TChunkList.Create;
  try

    load(ParamStr(1));

    makeChunks;

    kMeansReduce;

    dstData := srcData;

    save(ParamStr(2));

  finally
    chunkList.Free;
    Freemem(chunkAtts);
    Freemem(chunkFreqs);
    Freemem(srcData);
  end;
end.

