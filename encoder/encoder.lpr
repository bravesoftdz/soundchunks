program encoder;

{$mode objfpc}{$H+}

uses windows, Classes, sysutils, strutils, Types, fgl, MTProcs, math, extern, ap, conv, correlation, anysort;

const
  BandCount = 1;
  C1Freq = 32.703125;
  MaxChunksPerFrame = 4096;
  FrameLength = 10000; // im ms. if changed, adjust CLZRatio in PrepareFrames
  StreamVersion = 0;
  MaxAttenuation = 6;

type
  TEncoder = class;
  TFrame = class;
  TBand = class;
  TChunk = class;

  TBandGlobalData = record
    fcl, fch: Double;
    underSample: Integer;
    filteredData: TDoubleDynArray2;
    dstData: TSmallIntDynArray2;
  end;

  PBandGlobalData = ^TBandGlobalData;

  { TChunk }

  TChunkList = specialize TFPGObjectList<TChunk>;

  TChunk = class
  public
    frame: TFrame;
    reducedChunk: TChunk;

    channel, index, bandIndex: Integer;
    underSample: Integer;
    dstNegative: Boolean;
    dstAttenuation: Integer;

    origSrcData: PDouble;
    srcData: TDoubleDynArray;
    dct: TDoubleDynArray;
    dstData: TSmallIntDynArray;

    constructor Create(frm: TFrame; idx, bandIdx: Integer; underSmp: Integer; srcDta: PDouble);
    destructor Destroy; override;

    procedure ComputeDCT;
    procedure ComputeAttenuation;
    procedure MakeSrcData(origData: PDouble);
    procedure MakeDstData;
  end;

  { TBand }

  TBand = class
  public
    frame: TFrame;

    index: Integer;
    ChunkCount: Integer;

    srcData: array of PDouble;
    dstData: TDoubleDynArray2;

    finalChunks: TChunkList;

    globalData: PBandGlobalData;

    constructor Create(frm: TFrame; idx: Integer; startSample, endSample: Integer);
    destructor Destroy; override;

    procedure MakeChunks;
    procedure MakeDstData;
  end;

  { TFrame }

  TFrame = class
  public
    encoder: TEncoder;

    index: Integer;
    SampleCount: Integer;
    FrameSize: Integer;

    chunkRefs, reducedChunks: TChunkList;

    bands: array[0..BandCount - 1] of TBand;

    constructor Create(enc: TEncoder; idx, startSample, endSample: Integer);
    destructor Destroy; override;

    procedure MakeChunks;
    procedure KMeansReduce;
    procedure SaveStream(AStream: TStream);
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
    ChunkBitDepth: Integer; // 8 or 12 Bits
    ChunkSize: Integer;
    ChunksPerFrame: Integer;
    ReduceBassBand: Boolean;
    VariableFrameSizeRatio: Double;
    TrebleBoost: Boolean;
    ChunkBlend: Integer;

    ChannelCount: Integer;
    SampleRate: Integer;
    SampleCount: Integer;
    BlockSampleCount: Integer;
    ProjectedByteSize, FrameCount: Integer;
    Verbose: Boolean;

    srcHeader: array[$00..$2b] of Byte;
    srcData: TSmallIntDynArray2;
    dstData: TSmallIntDynArray2;

    frames: TFrameList;

    bandData: array[0 .. BandCount - 1] of TBandGlobalData;

    class function make16BitSample(smp: Double): SmallInt;
    class function makeOutputSample(smp: Double; OutBitDepth, Attenuation: Integer; Negative: Boolean): SmallInt;
    class function makeFloatSample(smp: SmallInt): Double; overload;
    class function makeFloatSample(smp: SmallInt; OutBitDepth, Attenuation: Integer; Negative: Boolean): Double; overload;
    class function ComputeAttenuation(chunkSz: Integer; const samples: TDoubleDynArray): Integer;
    class function ComputeDCT(chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
    class function ComputeInvDCT(chunkSz: Integer; const dct: TDoubleDynArray): TDoubleDynArray;
    class function ComputeDCT4(chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
    class function ComputeModifiedDCT(samplesSize: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
    class function ComputeInvModifiedDCT(dctSize: Integer; const dct: TDoubleDynArray): TDoubleDynArray;
    class function CompareEuclidean(firstCoeff, lastCoeff: Integer; const dctA, dctB: TDoubleDynArray): Double; overload;
    class function CompareEuclidean(firstCoeff, lastCoeff: Integer; const dctA, dctB: TSmallIntDynArray): Double; overload;
    class function CheckJoinPenalty(x, y, z, a, b, c: Double; TestRange: Boolean): Boolean; inline;
    class function ComputePsyADelta(const smpRef, smpTst: TSmallIntDynArray2): Double;
    class procedure createWAV(channels: word; resolution: word; rate: longint; fn: string; const data: TSmallIntDynArray);

    constructor Create(InFN, OutFN: String);
    destructor Destroy; override;

    procedure Load;
    procedure SaveWAV;
    procedure SaveGSC;
    procedure SaveStream(AStream: TStream);
    procedure SaveBandWAV(index: Integer; fn: String);

    procedure MakeBandGlobalData;
    procedure MakeBandSrcData(AIndex: Integer);

    procedure PrepareFrames;
    procedure MakeFrames;
    procedure MakeDstData;

    function DoFilterCoeffs(fc, transFactor: Double; HighPass, Windowed: Boolean): TDoubleDynArray;
    function DoFilter(const samples, coeffs: TDoubleDynArray): TDoubleDynArray;
    function DoBPFilter(fcl, fch, transFactor: Double; const samples: TDoubleDynArray): TDoubleDynArray;

    function ComputeEAQUAL(chunkSz: Integer; UseDIX, Verbz: Boolean; const smpRef, smpTst: TSmallIntDynArray): Double;
  end;


function IsDebuggerPresent(): LongBool stdcall; external 'kernel32.dll';

function HasParam(p: String): Boolean;
var i: Integer;
begin
  Result := False;
  for i := 3 to ParamCount do
    if SameText(p, ParamStr(i)) then
      Exit(True);
end;

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

function lerp(x, y, alpha: Double): Double; inline;
begin
  Result := x + (y - x) * alpha;
end;

function ilerp(x, y, alpha, maxAlpha: Integer): Integer; inline;
begin
  Result := x + ((y - x) * alpha) div maxAlpha;
end;

function revlerp(x, y, alpha: Double): Double; inline;
begin
  Result := (alpha - x) / (y - x);
end;

function nan0(x: Double): Double; inline;
begin
  Result := 0;
  if not IsNan(x) then
    Result := x;
end;

{ TChunk }

constructor TChunk.Create(frm: TFrame; idx, bandIdx: Integer; underSmp: Integer; srcDta: PDouble);
var
  i, j: Integer;
  mx: Double;
begin
  index := idx;
  bandIndex := bandIdx;
  underSample := underSmp;
  frame := frm;
  reducedChunk := Self;
  channel := -1;

  SetLength(srcData, frame.encoder.chunkSize);

  if Assigned(srcDta) then
  begin
    origSrcData := @srcDta[idx * (frame.encoder.chunkSize - frame.encoder.ChunkBlend) * underSample];
    MakeSrcData(origSrcData);
  end;

  // compute overall sign

  mx := 0.0;
  j := 0;
  for i := 0 to High(srcData) do
    if abs(srcData[i]) > mx then
    begin
      j := i;
      mx := abs(srcData[i]);
    end;
  dstNegative := (srcData[j] < 0);
end;

destructor TChunk.Destroy;
begin
  inherited Destroy;
end;

procedure TChunk.ComputeDCT;
var
  i: Integer;
  data: TDoubleDynArray;
begin
  SetLength(data, Length(dstData));
  for i := 0 to High(data) do
    data[i] := dstData[i] / ((1 shl (frame.encoder.ChunkBitDepth - 1)) - 1);
  dct := TEncoder.ComputeDCT4(Length(data), data);
end;

procedure TChunk.ComputeAttenuation;
begin
  dstAttenuation := TEncoder.ComputeAttenuation(Length(srcData), srcData);
end;

procedure TChunk.MakeSrcData(origData: PDouble);
var
  i, j, pos, n: Integer;
  f, acc: Double;
begin
  for i := 0 to High(srcData) do
  begin
    pos := i * underSample;

    acc := 0.0;
    n := 0;
    for j := 0 to underSample - 1 do
    begin
      if pos + j >= frame.SampleCount then
        Break;
      acc += origData[pos + j];
      Inc(n);
    end;

    if n = 0 then
      srcData[i] := 0
    else
      srcData[i] := acc / n;
  end;

  for i := 0 to frame.encoder.ChunkBlend - 1 do
  begin
    f := (i + 1) / (frame.encoder.ChunkBlend + 1);
    srcData[i] *= f;
    srcData[frame.encoder.ChunkSize - 1 - i] *= f;
  end;
end;

procedure TChunk.MakeDstData;
var
  i: Integer;
begin
  SetLength(dstData, length(srcData));
  for i := 0 to High(dstData) do
    dstData[i] := TEncoder.makeOutputSample(srcData[i], frame.encoder.ChunkBitDepth, dstAttenuation, dstNegative);
end;

{ TBand }

constructor TBand.Create(frm: TFrame; idx: Integer; startSample, endSample: Integer);
var
  i: Integer;
begin
  frame := frm;
  index := idx;
  globalData := @frame.encoder.bandData[index];

  SetLength(srcData, frame.encoder.ChannelCount);
  for i := 0 to High(srcData) do
    srcData[i] := @globalData^.filteredData[i, startSample];

  ChunkCount := (endSample - startSample + 1 - 1) div ((frame.encoder.ChunkSize - frame.encoder.ChunkBlend) * globalData^.underSample) + 1;

  finalChunks := TChunkList.Create;
end;

destructor TBand.Destroy;
begin
  finalChunks.Free;

  inherited Destroy;
end;

procedure TBand.MakeChunks;
var
  j, i: Integer;
  chunk: TChunk;
begin
  finalChunks.Clear;
  finalChunks.Capacity := ChunkCount * frame.encoder.ChannelCount;

  for i := 0 to ChunkCount - 1 do
    for j := 0 to frame.encoder.ChannelCount - 1 do
    begin
      chunk := TChunk.Create(frame, i, index, globalData^.underSample, srcData[j]);
      chunk.channel := j;
      chunk.ComputeAttenuation;
      chunk.MakeDstData;
      chunk.ComputeDCT;
      finalChunks.Add(chunk);
    end;
end;

procedure TBand.MakeDstData;
var
  i, j, k: Integer;
  chunk: TChunk;
  smp: Double;
  pos: TIntegerDynArray;
begin
  //WriteLn('MakeDstData #', index);

  SetLength(pos, frame.encoder.ChannelCount);
  SetLength(dstData, frame.encoder.ChannelCount, frame.SampleCount);
  for i := 0 to frame.encoder.ChannelCount - 1 do
  begin
    FillQWord(dstData[i, 0], frame.SampleCount, 0);
    pos[i] := 0;
  end;

  for i := 0 to finalChunks.Count - 1 do
  begin
    chunk := finalChunks[i];

    for j := 0 to frame.encoder.chunkSize - 1 do
    begin
      smp := TEncoder.makeFloatSample(chunk.reducedChunk.dstData[j], frame.encoder.ChunkBitDepth, chunk.dstAttenuation, chunk.dstNegative);

      for k := 0 to globalData^.underSample - 1 do
      begin
        if InRange(pos[chunk.channel], 0, High(dstData[chunk.channel])) then
          dstData[chunk.channel, pos[chunk.channel]] += smp;
        Inc(pos[chunk.channel]);
      end;
    end;

    Dec(pos[chunk.channel], frame.encoder.ChunkBlend);
  end;
end;

constructor TFrame.Create(enc: TEncoder; idx, startSample, endSample: Integer);
var
  i: Integer;
begin
  encoder := enc;
  index := idx;
  SampleCount := endSample - startSample + 1;

  for i := 0 to BandCount - 1 do
    bands[i] := TBand.Create(Self, i, startSample, endSample);

  if encoder.Verbose then
  begin
    Write('Frame #', index);
    for i := 0 to BandCount - 1 do
      Write(#9, bands[i].ChunkCount);
    WriteLn;
  end;

  chunkRefs := TChunkList.Create(False);
  reducedChunks := TChunkList.Create;
end;

destructor TFrame.Destroy;
var
  i: Integer;
begin
  reducedChunks.Free;
  chunkRefs.Free;

  for i := 0 to BandCount - 1 do
    bands[i].Free;

  inherited Destroy;
end;

procedure TFrame.MakeChunks;
var
  i, j, k: Integer;
begin
  chunkRefs.Clear;
  for i := Ord(not encoder.ReduceBassBand) to BandCount - 1 do
  begin
    bands[i].MakeChunks;
    for j := 0 to bands[i].finalChunks.Count - 1 do
      //for k := 1 to round(Power(bands[i].globalData^.underSample, sqrt(2.0))) do
        chunkRefs.Add(bands[i].finalChunks[j]);
  end;
end;

type
  TCountIndex = class
    Index, Count: Integer;
    Value: Double;
  end;

  TCountIndexList = specialize TFPGObjectList<TCountIndex>;

function CompareCountInvIndex (const Item1, Item2: TCountIndex): Integer;
begin
  Result := CompareValue(Item2.Count, Item1.Count);
end;

function CompareValueInvIndex (const Item1, Item2: TCountIndex): Integer;
begin
  Result := CompareValue(Item2.Value, Item1.Value);
end;

procedure TFrame.KMeansReduce;
var
  i, j, k, prec, colCount, clusterCount: Integer;
  chunk: TChunk;
  centroid: TDoubleDynArray;
  Clusters: TIntegerDynArray;
  Dataset: TFloatDynArray2;
  Centroids: TFloatDynArray2;
  Yakmo: PYakmo;
  CIList: TCountIndexList;
  CIInv: TIntegerDynArray;

  best: TANNFloat;
  bestIdx: Integer;
  tmp: TANNFloatDynArray2;
  tmp2: TANNFloatDynArray;
  KDT: PANNkdtree;
begin
  prec := encoder.Precision;

  colCount := encoder.chunkSize;
  clusterCount := encoder.ChunksPerFrame;

  SetLength(Dataset, chunkRefs.Count, colCount);

  for i := 0 to chunkRefs.Count - 1 do
    for j := 0 to colCount - 1 do
      Dataset[i, j] := chunkRefs[i].dct[j];

  if (prec > 0) and (chunkRefs.Count > clusterCount) then
  begin
    // usual chunk reduction using K-Means

    if encoder.Verbose then
      WriteLn('KMeansReduce Frame = ', index, ', N = ', chunkRefs.Count, ', K = ', clusterCount);

    SetLength(Clusters, chunkRefs.Count);
    SetLength(Centroids, clusterCount, colCount);

    Yakmo := yakmo_create(clusterCount, prec, MaxInt, 1, 0, 0, IfThen(encoder.Verbose, 1));
    yakmo_load_train_data(Yakmo, chunkRefs.Count, colCount, @Dataset[0]);
    yakmo_train_on_data(Yakmo, @Clusters[0]);
    yakmo_get_centroids(Yakmo, @Centroids[0]);
    yakmo_destroy(Yakmo);


    CIList := TCountIndexList.Create;
    try
      for i := 0 to clusterCount - 1 do
      begin
        CIList.Add(TCountIndex.Create);
        CIList[i].Index := i;

        for j := 0 to High(Clusters) do
          if Clusters[j] = i then
            Inc(CIList[i].Count);
      end;
      CIList.Sort(@CompareCountInvIndex);
      SetLength(CIInv, clusterCount);

      SetLength(centroid, encoder.chunkSize);

      reducedChunks.Clear;
      reducedChunks.Capacity := encoder.ChunksPerFrame;
      for i := 0 to clusterCount - 1 do
      begin
        chunk := TChunk.Create(Self, i, -1, 1, nil);
        reducedChunks.Add(chunk);

        for j := 0 to encoder.chunkSize - 1 do
          centroid[j] := nan0(Centroids[CIList[i].Index][j]);

        CIInv[CIList[i].Index] := i;

        centroid := TEncoder.ComputeDCT4(encoder.chunkSize, centroid);

        SetLength(chunk.dstData, encoder.chunkSize);
        for j := 0 to encoder.chunkSize - 1 do
          chunk.dstData[j] := EnsureRange(round(centroid[j] * ((1 shl (encoder.ChunkBitDepth - 1)) - 1)), -(1 shl (encoder.ChunkBitDepth - 1)), (1 shl (encoder.ChunkBitDepth - 1)) - 1);
  	  end;

      SetLength(tmp, clusterCount * MaxAttenuation * 2 {Negative}, encoder.chunkSize);
      SetLength(tmp2, encoder.chunkSize);
      for j := 0 to high(tmp) do
        for k := 0 to encoder.ChunkSize - 1 do
          tmp[j, k] := TEncoder.makeFloatSample(reducedChunks[(j shr 1) div MaxAttenuation].dstData[k], encoder.ChunkBitDepth, (j shr 1) mod MaxAttenuation, j and 1 <> 0);

      KDT := ann_kdtree_create(@tmp[0], Length(tmp), encoder.ChunkSize, 1, ANN_KD_SUGGEST);

      for i := 0 to chunkRefs.Count - 1 do
      begin
{$if true}
        for k := 0 to encoder.ChunkSize - 1 do
          tmp2[k] := chunkRefs[i].srcData[k];

        bestIdx := ann_kdtree_search(KDT, @tmp2[0], 0.0, @best);

        chunkRefs[i].dstNegative := bestIdx and 1 <> 0;
        chunkRefs[i].dstAttenuation := (bestIdx shr 1) mod MaxAttenuation;
        chunkRefs[i].reducedChunk := reducedChunks[(bestIdx shr 1) div MaxAttenuation];
{$else}
        chunkRefs[i].reducedChunk := reducedChunks[CIInv[Clusters[i]]];
{$ifend}
      end;

      ann_kdtree_destroy(KDT);
    finally
      CIList.Free;
    end;
  end
  else
  begin
    // passthrough mode

    reducedChunks.Clear;
    reducedChunks.Capacity := chunkRefs.Count;
    for i := 0 to reducedChunks.Capacity - 1 do
    begin
      chunk := TChunk.Create(Self, i, -1, 1, nil);

      reducedChunks.Add(chunk);

      chunk.dstData := Copy(chunkRefs[i].dstData);
    end;

    Centroids := Dataset;

    for i := 0 to chunkRefs.Count - 1 do
      chunkRefs[i].reducedChunk := reducedChunks[i];
  end;
end;

procedure TFrame.SaveStream(AStream: TStream);
var
  j, k, l, s1, s2: Integer;
  w : Integer;
  cl: TChunkList;
begin
  Assert(reducedChunks.Count <= MaxChunksPerFrame);

  w := (encoder.ChannelCount shl 8) or StreamVersion;
  AStream.WriteWord(w and $ffff);
  w := reducedChunks.Count or ((BandCount - 1) shl 13);
  AStream.WriteWord(w and $ffff);
  w := (encoder.ChunkSize shl 8) or encoder.ChunkBitDepth;
  AStream.WriteWord(w and $ffff);
  AStream.WriteDWord(encoder.SampleRate);

  cl := reducedChunks;
  if cl.Count = 0 then
    cl := chunkRefs;

  case encoder.ChunkBitDepth of
    8:
      for k := 0 to cl.Count - 1 do
        for l := 0 to encoder.ChunkSize - 1 do
          AStream.WriteByte((cl[k].dstData[l] - Low(ShortInt)) and $ff);
    12:
      for k := 0 to cl.Count - 1 do
        for l := 0 to encoder.ChunkSize div 2 - 1 do
        begin
          s1 := cl[k].dstData[l * 2 + 0] + 2048;
          s2 := cl[k].dstData[l * 2 + 1] + 2048;

          AStream.WriteByte(((s1 shr 4) and $f0) or ((s2 shr 8) and $0f));
          AStream.WriteByte(s1 and $ff);
          AStream.WriteByte(s2 and $ff);
        end;
    else
      Assert(False, 'ChunkBitDepth not supported');
  end;

  for j := 0 to BandCount - 1 do
  begin
    cl := bands[j].finalChunks;

    AStream.WriteDWord(cl.Count div encoder.ChannelCount);

    for k := 0 to cl.Count - 1 do
    begin
      w := cl[k].reducedChunk.index;
      w := w or (cl[k].dstAttenuation shl 12);
      w := w or IfThen(cl[k].dstNegative, $8000);
      AStream.WriteWord(w and $ffff);
    end;
  end;
end;

{ TEncoder }

procedure TEncoder.Load;
var
  wavFN: String;
  fs: TFileStream;
  i, j: Integer;
  data: TSmallIntDynArray;
begin
  if LowerCase(ExtractFileExt(inputFN)) <> '.wav' then
  begin
    WriteLn('Convert ', inputFN);
    wavFN := GetTempFileName + '.wav';
    DoExternalSOX(inputFN, wavFN);
  end
  else
  begin
    wavFN := inputFN;
  end;

  WriteLn('Load ', wavFN);

  fs := TFileStream.Create(wavFN, fmOpenRead or fmShareDenyNone);
  try
    fs.ReadBuffer(srcHeader[0], SizeOf(srcHeader));
    SampleRate := PInteger(@srcHeader[$18])^;
    ChannelCount := PWORD(@srcHeader[$16])^;

    SampleCount := (fs.Size - fs.Position) div (SizeOf(SmallInt) * ChannelCount);
    SetLength(srcData, ChannelCount, SampleCount);

    SetLength(data, SampleCount * ChannelCount);
    fs.ReadBuffer(data[0], SampleCount * ChannelCount * 2);

    for i := 0 to SampleCount - 1 do
      for j := 0 to ChannelCount - 1 do
        srcData[j, i] := data[i * ChannelCount + j];
  finally
    fs.Free;

    if wavFN <> inputFN then
      DeleteFile(wavFN);
  end;
end;

procedure TEncoder.SaveWAV;
var
  i, j: Integer;
  fs: TFileStream;
  wavFN: String;
  data: TSmallIntDynArray;
begin
  wavFN := ChangeFileExt(outputFN, '.wav');

  WriteLn('Save ', wavFN);

  fs := TFileStream.Create(wavFN, fmCreate or fmShareDenyWrite);
  try
    fs.WriteBuffer(srcHeader[0], SizeOf(srcHeader));

    SetLength(data, SampleCount * ChannelCount);

    for i := 0 to SampleCount - 1 do
      for j := 0 to ChannelCount - 1 do
        data[i * ChannelCount + j] := dstData[j, i];

    fs.WriteBuffer(data[0], SampleCount * ChannelCount * 2);
  finally
    fs.Free;
  end;
end;

procedure TEncoder.SaveGSC;
var
  fs: TFileStream;
  cur: TMemoryStream;
  fn: String;
begin
  fs := nil;
  fn := ChangeFileExt(outputFN, '.gsc');
  cur := TMemoryStream.Create;
  fs := TFileStream.Create(fn, fmCreate or fmShareDenyWrite);
  try
    WriteLn('Save ', fn);

    SaveStream(cur);
    cur.Position := 0;

    fs.CopyFrom(cur, cur.Size);

    writeln('FinalByteSize = ', cur.Size);
    writeln('FinalBitRate = ', round(cur.size * (8 / 1024) / (SampleCount / SampleRate)));
  finally
    fs.Free;
    cur.Free;
  end;
end;

procedure TEncoder.SaveStream(AStream: TStream);
var
  i: Integer;
  ZStream: TMemoryStream;
begin
  ZStream := TMemoryStream.Create;
  try
    for i := 0 to FrameCount - 1 do
      frames[i].SaveStream(ZStream);

    LZCompress(ZStream, False, False, AStream);
  finally
    ZStream.Free;
  end;
end;

procedure TEncoder.SaveBandWAV(index: Integer; fn: String);
var
  i, j: Integer;
  fs: TFileStream;
  data: TSmallIntDynArray;
begin
  //WriteLn('SaveBandWAV #', index, ' ', fn);

  fs := TFileStream.Create(fn, fmCreate or fmShareDenyWrite);
  try
    fs.WriteBuffer(srcHeader[0], SizeOf(srcHeader));

    SetLength(data, SampleCount * ChannelCount);

    for i := 0 to SampleCount - 1 do
      for j := 0 to ChannelCount - 1 do
        data[i * ChannelCount + j] := bandData[index].dstData[j, i];

    fs.WriteBuffer(data[0], SampleCount * ChannelCount * 2);
  finally
    fs.Free;
  end;
end;

procedure TEncoder.MakeBandGlobalData;
var
  i: Integer;
  ratio, hc: Double;
  bnd: TBandGlobalData;
begin
  for i := 0 to BandCount - 1 do
  begin
    bnd.dstData := nil;
    FillChar(bnd, SizeOf(bnd), 0);

    // determing low and high bandpass frequencies

    hc := min(HighCut, SampleRate / 2);
    ratio := (log2(hc) - log2(max(C1Freq, LowCut))) / BandCount;

    if i = 0 then
      bnd.fcl := LowCut / SampleRate
    else
      bnd.fcl := 0.5 * power(2.0, -floor((BandCount - i) * ratio));

    if i = BandCount - 1 then
      bnd.fch := hc / SampleRate
    else
      bnd.fch := 0.5 * power(2.0, -floor((BandCount - 1 - i) * ratio));

    // undersample if the band high freq is a lot lower than nyquist

    bnd.underSample := Max(1, round(0.25 / bnd.fch));

    bandData[i] := bnd;
  end;
end;

procedure TEncoder.MakeBandSrcData(AIndex: Integer);
var
  i, j: Integer;
  bnd: TBandGlobalData;
begin
  bnd := bandData[AIndex];

  SetLength(bnd.filteredData, ChannelCount, SampleCount);
  for i := 0 to ChannelCount - 1 do
    for j := 0 to SampleCount - 1 do
      bnd.filteredData[i, j] := makeFloatSample(srcData[i, j]);

  // band pass the samples
  for i := 0 to ChannelCount - 1 do
    bnd.filteredData[i] := DoBPFilter(bnd.fcl, bnd.fch, BandTransFactor, bnd.filteredData[i]);

  bandData[AIndex] := bnd;
end;

procedure TEncoder.PrepareFrames;

  procedure DoBand(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  begin
    MakeBandSrcData(AIndex);
  end;

const
  CLZRatio = 0.92;
var
  j, i, k, nextStart, psc, tentativeByteSize: Integer;
  frm: TFrame;
  fixedCost, frameCost, bandCost, avgPower, totalPower, perFramePower, curPower, smp: Double;
begin
  MakeBandGlobalData;

  // pass 1

  BlockSampleCount := 0;
  for i := 0 to BandCount - 1 do
    if bandData[i].underSample * (ChunkSize - ChunkBlend) > BlockSampleCount then
      BlockSampleCount := bandData[i].underSample * (ChunkSize - ChunkBlend);

  // ensure srcData ends on a full block
  psc := SampleCount;
  SampleCount := ((SampleCount - 1) div BlockSampleCount + 1) * BlockSampleCount;
  SetLength(srcData, ChannelCount, SampleCount);
  for j := 0 to ChannelCount - 1 do
    for i := psc to SampleCount - 1 do
      srcData[j, i] := 0;

  if BitRate > 0 then
    ProjectedByteSize := ceil((SampleCount / SampleRate) * (BitRate * 1024 / 8))
  else
    ProjectedByteSize := MaxInt;

  if Verbose then
  begin
    writeln('ProjectedByteSize = ', ProjectedByteSize);
  end;

  FrameCount := ceil(SampleCount / (SampleRate * (FrameLength / 1000)));

  Inc(ChunksPerFrame);
  repeat
    Dec(ChunksPerFrame);

    fixedCost := 0 {no header besides frame};

    bandCost := 0;
    for i := 0 to BandCount - 1 do
      bandCost += (SampleCount * ChannelCount * (Log2(ChunksPerFrame) + Log2(MaxAttenuation) + 1 {dstNegative})) / (8 {bytes -> bits} * (ChunkSize - ChunkBlend) * bandData[i].underSample);

    frameCost := (ChunksPerFrame * ChunkSize) * ChunkBitDepth / 8 + (3 * SizeOf(Word) + (1 + BandCount) * SizeOf(Cardinal)) {frame header};

    tentativeByteSize := Round((fixedCost + bandCost + FrameCount * frameCost) * CLZRatio);

  until (tentativeByteSize <= ProjectedByteSize) or (ChunksPerFrame <= 1);

  ProjectedByteSize := tentativeByteSize;

  writeln('ChannelCount = ', ChannelCount);
  writeln('SampleRate = ', SampleRate);
  writeln('FrameCount = ', FrameCount);
  writeln('ChunksPerFrame = ', ChunksPerFrame);

  Assert(ChunksPerFrame > 0, 'Null ChunksPerFrame! (BitRate too low)');

  if Verbose then
  begin
    WriteLn('SampleCount = ', SampleCount);
    writeln('FrameSize = ', ProjectedByteSize div FrameCount);
    writeln('ProjectedByteSize = ', ProjectedByteSize);
    writeln('ChunkSize = ', ChunkSize);
  end;

  ProcThreadPool.DoParallelLocalProc(@DoBand, 0, BandCount - 1, nil);

  // pass 2

  avgPower := 0.0;
  for j := 0 to ChannelCount - 1 do
    for i := 0 to SampleCount - 1 do
      avgPower += Sqr(makeFloatSample(srcData[j, i]));
  avgPower := Sqrt(avgPower) / (SampleCount * ChannelCount);

  totalPower := 0.0;
  for i := 0 to SampleCount - 1 do
  begin
    smp := 0.0;
    for j := 0 to ChannelCount - 1 do
      smp += Sqr(makeFloatSample(srcData[j, i]));
    smp := Sqrt(smp) / ChannelCount;

    totalPower += 1.0 - lerp(avgPower, smp, VariableFrameSizeRatio);
  end;

  perFramePower := totalPower / FrameCount;

  if Verbose then
  begin
    writeln('TotalPower = ', FormatFloat('0.00', totalPower));
    writeln('PerFramePower = ', FormatFloat('0.00', perFramePower));
  end;

  k := 0;
  nextStart := 0;
  curPower := 0.0;
  for i := 0 to SampleCount - 1 do
  begin
    smp := 0.0;
    for j := 0 to ChannelCount - 1 do
      smp += Sqr(makeFloatSample(srcData[j, i]));
    smp := Sqrt(smp) / ChannelCount;

    curPower += 1.0 - lerp(avgPower, smp, VariableFrameSizeRatio);

    if (i mod BlockSampleCount = 0) and (curPower >= perFramePower) then
    begin
      frm := TFrame.Create(Self, k, nextStart, i - 1);
      frames.Add(frm);

      curPower := 0.0;
      nextStart := i;
      Inc(k);
    end;
  end;

  frm := TFrame.Create(Self, k, nextStart, SampleCount - 1);
  frames.Add(frm);

  FrameCount := frames.Count;

  for i := 0 to BandCount - 1 do
     WriteLn('Band #', i, ' (', round(bandData[i].fcl * SampleRate), ' Hz .. ', round(bandData[i].fch * SampleRate), ' Hz); ', bandData[i].underSample);
end;

procedure TEncoder.MakeFrames;

  procedure DoFrame(AIndex: PtrInt; AData: Pointer; AItem: TMultiThreadProcItem);
  var
    i: Integer;
    frm: TFrame;
  begin
    frm := frames[AIndex];

    frm.MakeChunks;
    frm.KMeansReduce;
    for i := 0 to BandCount - 1 do
      frm.bands[i].MakeDstData;
    Write('.');
  end;
begin
  ProcThreadPool.DoParallelLocalProc(@DoFrame, 0, FrameCount - 1, nil);
  WriteLn;
end;

function TEncoder.DoFilter(const samples, coeffs: TDoubleDynArray): TDoubleDynArray;
var
  i: Integer;
begin
  Result := nil;
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
    coeffs := DoFilterCoeffs(fcl, transFactor * fcl, True, True);
    Result := DoFilter(Result, coeffs);
  end;

  if fch < 0.5 then
  begin
    coeffs := DoFilterCoeffs(fch, transFactor * fch, False, True);
    Result := DoFilter(Result, coeffs);
  end;
end;

constructor TEncoder.Create(InFN, OutFN: String);
begin
  inputFN := InFN;
  outputFN := OutFN;

  BitRate := -1;
  Precision := 1;
  LowCut := 0.0;
  HighCut := 24000.0;
  ChunkBitDepth := 12;
  ChunkSize := 4;
  ReduceBassBand := True;
  TrebleBoost := False;
  VariableFrameSizeRatio := 0.0;
  ChunkBlend := 0;

  ChunksPerFrame := MaxChunksPerFrame;
  BandTransFactor := 1 / 256;

  frames := TFrameList.Create;
end;

destructor TEncoder.Destroy;
begin
  frames.Free;

  inherited Destroy;
end;

procedure TEncoder.MakeDstData;
var
  i, j, k, l, pos: Integer;
  smp: Double;
  bnd: TBandGlobalData;
  resamp: array[0 .. BandCount - 1] of TDoubleDynArray;
  floatDst: TDoubleDynArray;
begin
  WriteLn('MakeDstData');

  SetLength(dstData, ChannelCount, SampleCount);
  for l := 0 to ChannelCount - 1 do
    FillWord(dstData[l, 0], Length(dstData), 0);

  SetLength(floatDst, SampleCount);

  for i := 0 to BandCount - 1 do
  begin
    bnd := bandData[i];

    SetLength(bnd.dstData, ChannelCount, SampleCount);
    for l := 0 to ChannelCount - 1 do
      FillWord(bnd.dstData[l, 0], SampleCount, 0);

    bandData[i] := bnd;
  end;

  for l := 0 to ChannelCount - 1 do
  begin
    FillQWord(floatDst[0], Length(floatDst), 0);

    pos := 0;
    for k := 0 to frames.Count - 1 do
    begin
      for j := 0 to BandCount - 1 do
      begin
        bnd := bandData[j];
        resamp[j] := DoBPFilter(bnd.fcl, bnd.fch, BandTransFactor, frames[k].bands[j].dstData[l]);
      end;

      for i := 0 to frames[k].SampleCount - 1 do
      begin
        for j := 0 to BandCount - 1 do
        begin
{$if true}
          smp := frames[k].bands[j].dstData[l, i];
{$else}
          smp := resamp[j][i];
{$endif}

          if InRange(pos, 0, High(dstData[l])) then
          begin
            bandData[j].dstData[l, pos] := make16BitSample(smp);
            floatDst[pos] := floatDst[pos] + smp;
          end;
        end;

        Inc(pos);
      end;
    end;

    for i := 0 to High(floatDst) do
      dstData[l, i] := make16BitSample(floatDst[i]);
  end;
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
{$if true}
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
  Result := EnsureRange(round(smp * High(SmallInt)), Low(SmallInt), High(SmallInt));
end;

class function TEncoder.makeFloatSample(smp: SmallInt): Double;
begin
  Result := smp / High(SmallInt);
end;

class function TEncoder.makeOutputSample(smp: Double; OutBitDepth, Attenuation: Integer; Negative: Boolean): SmallInt;
var
  obd: Integer;
  smp16: SmallInt;
begin
  obd := (1 shl (OutBitDepth - 1));
  smp16 := round(smp * obd * (1 + Attenuation));
  if Negative then smp16 := -smp16;
  smp16 := EnsureRange(smp16, -obd + 1, obd - 1);
  Result := smp16;
end;

class function TEncoder.makeFloatSample(smp: SmallInt; OutBitDepth, Attenuation: Integer; Negative: Boolean): Double;
var
  obd: Integer;
  smp16: SmallInt;
begin
  obd := (1 shl (OutBitDepth - 1)) * (1 + Attenuation);
  smp16 := smp;
  if Negative then smp16 := -smp16;
  Result := smp16 / obd;
  Result := EnsureRange(Result, -1.0, 1.0);
end;

class function TEncoder.ComputeAttenuation(chunkSz: Integer; const samples: TDoubleDynArray): Integer;
var
  i, hiSmp: Integer;
begin
  hiSmp := 0;
  for i := 0 to chunkSz - 1 do
    hiSmp := max(hiSmp, ceil(abs(samples[i] * High(SmallInt))));

  Result := 1;
  repeat
    Inc(Result)
  until (hiSmp * Result > High(SmallInt)) or (Result > MaxAttenuation);
  Dec(Result, 2);
end;

class function TEncoder.ComputeDCT(chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
var
  k, n: Integer;
  sum, s: Double;
begin
  SetLength(Result, length(samples));
  for k := 0 to chunkSz - 1 do
  begin
    s := ifthen(k = 0, sqrt(0.5), 1.0);

    sum := 0;
    for n := 0 to chunkSz - 1 do
      sum += s * samples[n] * cos(pi / chunkSz * (n + 0.5) * k);

    Result[k] := sum * sqrt (2.0 / chunkSz);
  end;
end;

class function TEncoder.ComputeInvDCT(chunkSz: Integer; const dct: TDoubleDynArray): TDoubleDynArray;
var
  k, n: Integer;
  sum: Double;
begin
  SetLength(Result, length(dct));
  for k := 0 to chunkSz - 1 do
  begin
    sum := sqrt(0.5) * dct[0];
    for n := 1 to chunkSz - 1 do
      sum += dct[n] * cos (pi / chunkSz * (k + 0.5) * n);

    Result[k] := sum * sqrt (2.0 / chunkSz);
  end;
end;

class function TEncoder.ComputeDCT4(chunkSz: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
var
  k, n: Integer;
  sum: Double;
begin
  SetLength(Result, length(samples));
  for k := 0 to chunkSz - 1 do
  begin
    sum := 0;
    for n := 0 to chunkSz - 1 do
      sum += samples[n] * cos(pi / chunkSz * (n + 0.5) * (k + 0.5));

    Result[k] := sum * sqrt (2.0 / chunkSz);
  end;
end;

// MDCT cannot be used (would need overlapped add in decoder)
class function TEncoder.ComputeModifiedDCT(samplesSize: Integer; const samples: TDoubleDynArray): TDoubleDynArray;
var
  k, n: Integer;
  sum: Double;
begin
  SetLength(Result, length(samples) div 2);
  for k := 0 to samplesSize div 2 - 1 do
  begin
    sum := 0;
    for n := 0 to samplesSize - 1 do
      sum += samples[n] * cos(pi / (samplesSize div 2) * (n + 0.5 + (samplesSize div 2) * 0.5) * (k + 0.5));

    Result[k] := sum;
  end;
end;

// IMDCT cannot be used (would need overlapped add in decoder)
class function TEncoder.ComputeInvModifiedDCT(dctSize: Integer; const dct: TDoubleDynArray): TDoubleDynArray;
var
  k, n, i: Integer;
  sum: Double;
begin
  SetLength(Result, length(dct));
  for n := 0 to dctSize - 1 do
  begin
    sum := 0;
    for k := 0 to dctSize div 2 - 1 do
      sum += dct[k] * cos (pi / (dctSize div 2) * (n + 0.5 + (dctSize div 2) * 0.5) * (k + 0.5));

    Result[n] := sum / (dctSize div 2);
  end;

  for i := 0 to dctSize div 2 - 1 do
  begin
    Result[i] += Result[i + dctSize div 2];
    Result[i + dctSize div 2] := 0.0;
  end;
end;

class function TEncoder.CompareEuclidean(firstCoeff, lastCoeff: Integer; const dctA, dctB: TDoubleDynArray): Double;
var
  i: Integer;
begin
  Result := 0.0;

  for i := firstCoeff to lastCoeff do
    Result += sqr(dctA[i] - dctB[i]);

  Result := sqrt(Result) / (lastCoeff - firstCoeff + 1);
end;

class function TEncoder.CompareEuclidean(firstCoeff, lastCoeff: Integer; const dctA, dctB: TSmallIntDynArray): Double;
var
  i: Integer;
begin
  Result := 0.0;

  for i := firstCoeff to lastCoeff do
    Result += sqr((dctA[i] - dctB[i]) / High(SmallInt));

  Result := sqrt(Result) / (lastCoeff - firstCoeff + 1);
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
  FNTmp, FNRef, FNTst: String;
begin
  FNTmp := GetTempFileName('', 'tmp-'+IntToStr(GetCurrentThreadId))+'.wav';
  FNRef := GetTempFileName('', 'ref-'+IntToStr(GetCurrentThreadId))+'.wav';
  FNTst := GetTempFileName('', 'tst-'+IntToStr(GetCurrentThreadId))+'.wav';

  createWAV(ChannelCount, 16, SampleRate, FNTmp, smpRef);
  DoExternalSOX(FNTmp, FNRef, 48000);

  createWAV(ChannelCount, 16, SampleRate, FNTmp, smpTst);
  DoExternalSOX(FNTmp, FNTst, 48000);

  Result := DoExternalEAQUAL(FNRef, FNTst, Verbz, UseDIX, -1);

  DeleteFile(FNTst);
  DeleteFile(FNRef);
  DeleteFile(FNTst);
end;

class function TEncoder.ComputePsyADelta(const smpRef, smpTst: TSmallIntDynArray2): Double;
var
  i, j, len: Integer;
  rr, rt: TReal1DArray;
begin
  len := length(smpRef) * length(smpRef[0]);
  Assert(len = length(smpTst) * length(smpTst[0]), 'ComputePsyADelta length mismatch!');
  SetLength(rr, len);
  SetLength(rt, len);

  for j := 0 to High(smpRef) do
    for i := 0 to High(smpRef[0]) do
    begin
      rr[j * Length(smpRef[0]) + i] := smpRef[j, i];
      rt[j * Length(smpRef[0]) + i] := smpTst[j, i];
    end;

  //rr := ComputeDCT(len, rr);
  //rt := ComputeDCT(len, rt);

  Result := CompareEuclidean(0, len - 1, rr, rt) * length(smpRef);
end;

class procedure TEncoder.createWAV(channels: word; resolution: word; rate: longint; fn: string; const data: TSmallIntDynArray);
var
  wf : TFileStream;
  wh : TWavHeader;
begin
  wh.rId             := $46464952; { 'RIFF' }
  wh.rLen            := 36 + Length(data) * SizeOf(data[0]); { length of sample + format }
  wh.wId             := $45564157; { 'WAVE' }
  wh.fId             := $20746d66; { 'fmt ' }
  wh.fLen            := 16; { length of format chunk }
  wh.wFormatTag      := 1; { PCM data }
  wh.nChannels       := channels; { mono/stereo }
  wh.nSamplesPerSec  := rate; { sample rate }
  wh.nAvgBytesPerSec := channels*rate*(resolution div 8);
  wh.nBlockAlign     := channels*(resolution div 8);
  wh.wBitsPerSample  := resolution;{ resolution 8/16 }
  wh.dId             := $61746164; { 'data' }
  wh.wSampleLength   := Length(data) * SizeOf(data[0]); { sample size }

  wf := TFileStream.Create(fn, fmCreate or fmShareDenyNone);
  try
    wf.WriteBuffer(wh, SizeOf(wh));
    wf.WriteBuffer(data[0], Length(data) * SizeOf(data[0]));
  finally
    wf.Free;
  end;
end;


procedure test_makeSample;
var
  i: Integer;
  smp, o, so: SmallInt;
  obd, bs: SmallInt;
  sgn: Boolean;
  f, sf: Double;
begin
  for i := 0 to 65535 do
  begin
    bs := RandomRange(0, 7);
    sgn := Random >= 0.5;
    obd := 12;//RandomRange(1, 8);

    smp := (i mod (1 shl obd)) - (1 shl (obd - 1));
    sf := smp / (1 shl (obd - 1)) / (1 * (1 + bs)) * IfThen(sgn, -1, 1);

    f := TEncoder.makeFloatSample(smp, obd, bs, sgn);
    o := TEncoder.makeOutputSample(f, obd, bs, sgn);
    so := TEncoder.makeOutputSample(sf, obd, bs, sgn);
    writeln(smp,#9,o,#9,so,#9,bs,#9,sgn,#9,FloatToStr(f));
    assert(smp = o);
    assert(smp = so);
  end;

  halt;
end;

var
  enc: TEncoder;
  i: Integer;
  dix, psy: double;
  s: String;
begin
  try
    FormatSettings.DecimalSeparator := '.';

{$ifdef DEBUG}
    //ProcThreadPool.MaxThreadCount := 1;
{$else}
    SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS);
{$endif}

    //test_makeSample;


    if ParamCount < 2 then
    begin
      WriteLn('Usage: ', ExtractFileName(ParamStr(0)) + ' <source file> <dest file> [options]');
      Writeln('Main options:');
      WriteLn(#9'-br'#9'encoder bit rate in kilobits/second; example: "-br320"');
      WriteLn(#9'-lc'#9'bass cutoff frequency');
      WriteLn(#9'-hc'#9'treble cutoff frequency');
      WriteLn(#9'-vfr'#9'RMS power based variable frame size ratio (0.0-1.0); default: "-vfr1.0"');
      WriteLn(#9'-v'#9'verbose mode');
      Writeln('Development options:');
      WriteLn(#9'-cs'#9'chunk size');
      WriteLn(#9'-cpf'#9'max. chunks per frame (256-4096)');
      WriteLn(#9'-pbb'#9'disable lossy compression on bass band');
      WriteLn(#9'-cbd'#9'chunk bit depth (8,12)');
      WriteLn(#9'-pr'#9'K-means precision; 0: "lossless" mode');
      WriteLn(#9'-cb'#9'chunk blend');

      WriteLn;
      Writeln('(source file must be 16bit WAV or anything SOX can convert)');
      WriteLn;
      Exit;
    end;

    enc := TEncoder.Create(ParamStr(1), ParamStr(2));
    try
      enc.BitRate := round(ParamValue('-br', enc.BitRate));
      enc.Precision := round(ParamValue('-pr', enc.Precision));
      enc.LowCut := ParamValue('-lc', enc.LowCut);
      enc.HighCut := ParamValue('-hc', enc.HighCut);
      enc.VariableFrameSizeRatio :=  EnsureRange(ParamValue('-vfr', enc.VariableFrameSizeRatio), 0.0, 1.0);
      enc.ChunkBitDepth := EnsureRange(round(ParamValue('-cbd', enc.ChunkBitDepth)), 1, 16);
      enc.ChunkSize := round(ParamValue('-cs', enc.ChunkSize));
      enc.ChunksPerFrame := EnsureRange(round(ParamValue('-cpf', enc.ChunksPerFrame)), 256, MaxChunksPerFrame);
      enc.Verbose := HasParam('-v');
      enc.ReduceBassBand := not HasParam('-pbb');
      enc.ChunkBlend := EnsureRange(round(ParamValue('-cb', enc.ChunkBlend)), 0, enc.ChunkSize div 2);

      WriteLn('BitRate = ', FloatToStr(enc.BitRate));
      WriteLn('LowCut = ', FloatToStr(enc.LowCut));
      WriteLn('HighCut = ', FloatToStr(enc.HighCut));
      WriteLn('VariableFrameSizeRatio = ', FloatToStr(enc.VariableFrameSizeRatio));
      if enc.Verbose then
      begin
        WriteLn('ChunkSize = ', enc.ChunkSize);
        WriteLn('MaxChunksPerFrame = ', enc.ChunksPerFrame);
        WriteLn('ReduceBassBand = ', BoolToStr(enc.ReduceBassBand, True));
        WriteLn('ChunkBitDepth = ', enc.ChunkBitDepth);
        WriteLn('Precision = ', enc.Precision);
        WriteLn('ChunkBlend = ', enc.ChunkBlend);
      end;
      WriteLn;

      enc.Load;

      enc.PrepareFrames;
      enc.MakeFrames;
      enc.MakeDstData;

      enc.SaveWAV;
      if BandCount > 1 then
        for i := 0 to BandCount - 1 do
          enc.SaveBandWAV(i, ChangeFileExt(enc.outputFN, '-' + IntToStr(i) + '.wav'));
      if enc.Precision > 0 then
        enc.SaveGSC;

      //dix := enc.ComputeEAQUAL(enc.SampleCount, True, True, enc.srcData, enc.dstData);
      //WriteLn('EAQUAL = ', FloatToStr(dix));

      psy := enc.ComputePsyADelta(enc.srcData, enc.dstData);
      WriteLn('PsyADelta = ', FormatFloat(',0.0000000000', psy));

      s := FloatToStr(dix) + ' ' + FormatFloat(',0.0000000000', psy) + ' ';
      for i := 0 to ParamCount do s := s + ParamStr(i) + ' ';
      ShellExecute(0, 'open', 'cmd.exe', PChar('/c echo ' + s + ' >> ..\log.txt'), '', 0);

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

