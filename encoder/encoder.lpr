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
    reducedData: PSmallInt;
    dct: TSingleDynArray;
    weightedDct: TSingleDynArray;

    constructor Create(enc: TEncoder; idx: Integer);

    procedure ComputeDCT;
    procedure ComputeInvDCT;
  end;

  TChunkList = specialize TFPGObjectList<TChunk>;

  { TEncoder }

  TEncoder = class
  private
    procedure kMeansReduce;
    procedure load(fn: String);
    procedure makeChunks;
    procedure makeDstData;
    procedure save(fn: String);
  public
    chunkSize: Integer;
    quality: Single;
    desiredChunkCount: Integer;
    restartCount: Integer;
    sampleRate: Integer;
    joinSize: Integer;

    srcHeader: array[$00..$2b] of Byte;
    srcData: TSmallIntDynArray;
    srcDataCount: Integer;
    dstData: TSmallIntDynArray;
    chunkList: TChunkList;
    chunkFreqs: TSingleDynArray;
    chunkWgtAtts: TSingleDynArray;
    chunkJoinAtts: TSingleDynArray;

    constructor Create;
    destructor Destroy; override;
  end;


function make16BitSample(smp: Single): SmallInt;
begin
  Result := EnsureRange(round(smp), Low(SmallInt), High(SmallInt));
end;


{ TChunk }

constructor TChunk.Create(enc: TEncoder; idx: Integer);
begin
  encoder := enc;
  SetLength(rawData, encoder.chunkSize);
  Move(encoder.srcData[idx * (encoder.chunkSize - encoder.joinSize)], rawData[0], encoder.chunkSize * SizeOf(SmallInt));
  reducedData := @rawData[0];
  SetLength(dct, encoder.chunkSize);
  SetLength(weightedDct, encoder.chunkSize);
end;

procedure TChunk.ComputeDCT;
var
  k, n: Integer;
  sum, s: Single;
begin
  for k := 0 to encoder.chunkSize - 1 do
  begin
    sum := 0;
    s := ifthen(k = 0, sqrt(0.5), 1.0);
    for n := 0 to encoder.chunkSize - 1 do
      sum += s * rawData[n] * cos(pi * (n + 0.5) * k / encoder.chunkSize);
    dct[k] := sum * sqrt (2.0 / encoder.chunkSize);
    weightedDct[k] := dct[k] * encoder.chunkWgtAtts[k];
  end;
end;

procedure TChunk.ComputeInvDCT;
var
  k, n: Integer;
  sum, s: Single;
begin
  for n := 0 to encoder.chunkSize - 1 do
  begin
    sum := 0;
    for k := 0 to encoder.chunkSize - 1 do
    begin
      s := ifthen(k = 0, sqrt(0.5), 1.0);
      sum += s * dct[k] * cos (pi * (n + 0.5) * k / encoder.chunkSize);
    end;
    rawData[n] := make16BitSample(sum * sqrt(2.0 / encoder.chunkSize));
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
end;

destructor TEncoder.Destroy;
begin
  chunkList.Free;

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
  i: Integer;
  chunk: TChunk;
  chunkCount: Integer;
  fsq, att: Single;
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
    att := maxvalue([0.0, joinSize - i, i - chunkSize + joinSize]);
    if joinSize = 0 then
    begin
      chunkJoinAtts[i] := 1.0;
    end
    else
    begin
{$if false}
      chunkJoinAtts[i] := (joinSize - att) / joinSize;
{$else}
      chunkJoinAtts[i] := 0.5 - cos((joinSize - att) / joinSize * pi) * 0.5;
{$ifend}
    end;
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
        if XYC[i] = j then
          chunkList[i].reducedData := @chunkList[First].rawData[0];
  end;
end;

procedure TEncoder.makeDstData;
var
  i, j, pos: Integer;
  smp: Single;
begin
  WriteLn('makeDstData');

  SetLength(dstData, Length(srcData));
  FillWord(dstData[0], Length(srcData), 0);
  for i := 0 to chunkList.Count - 1 do
    for j := 0 to chunkSize - 1 do
    begin
      pos := i * (chunkSize - joinSize) + j;
      smp := dstData[pos] + chunkJoinAtts[j] * chunkList[i].reducedData[j];
      dstData[pos] := make16BitSample(smp);
    end;
end;


var
  enc: TEncoder;

begin
  FormatSettings.DecimalSeparator := '.';

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
end.

