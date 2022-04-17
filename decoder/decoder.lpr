program decoder;

uses Types, SysUtils, Classes, Math, extern;

const
  CAttrMul = round(High(SmallInt) * (High(SmallInt) / 2047));
  CAttrLookup : array[0 .. 15] of Integer = (
    CAttrMul div (1 + 0),
    CAttrMul div (1 + 1),
    CAttrMul div (1 + 2),
    CAttrMul div (1 + 3),
    CAttrMul div (1 + 4),
    CAttrMul div (1 + 5),
    CAttrMul div (1 + 6),
    CAttrMul div (1 + 7),
    -CAttrMul div (1 + 0),
    -CAttrMul div (1 + 1),
    -CAttrMul div (1 + 2),
    -CAttrMul div (1 + 3),
    -CAttrMul div (1 + 4),
    -CAttrMul div (1 + 5),
    -CAttrMul div (1 + 6),
    -CAttrMul div (1 + 7)
  );


  function CreateWAVHeader(channels: word; resolution: word; rate, size: longint): TWavHeader;
  var
    wh : TWavHeader;
  begin
    wh.rId             := $46464952; { 'RIFF' }
    wh.rLen            := 36 + size; { length of sample + format }
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
    wh.wSampleLength   := size; { sample size }

    Result := wh;
  end;

  procedure GSCUnpack(ASourceStream, ADestStream: TStream);
  var
    i, j, k, b, s1, s2: Integer;
    w: Word;
    chunkIndex, chunkAttrs: TIntegerDynArray;
    StreamVersion, ChannelCount, ChunkBitDepth, ChunkSize, ChunkCount, FrameLength, SampleRate: Integer;
    Chunks: TSmallIntDynArray2;
    memStream: TMemoryStream;
  begin
    memStream := TMemoryStream.Create;
    try
      while ASourceStream.Position <> ASourceStream.Size do
      begin
        StreamVersion := ASourceStream.ReadByte;
        ChannelCount := ASourceStream.ReadByte;
        ChunkCount := ASourceStream.ReadWord and $1fff;
        ChunkBitDepth := ASourceStream.ReadByte;
        ChunkSize := ASourceStream.ReadByte;
        SampleRate := ASourceStream.ReadDWord;

        SetLength(Chunks, ChunkCount, ChunkSize);

        case ChunkBitDepth of
          8:
            for i := 0 to ChunkCount - 1 do
              for j := 0 to ChunkSize - 1 do
              begin
                b := ASourceStream.ReadByte;
                Chunks[i, j] := (b + Low(ShortInt)) * 2047 div High(ShortInt);
              end;
          12:
            for i := 0 to ChunkCount - 1 do
              for j := 0 to ChunkSize div 2 - 1 do
              begin
                b := ASourceStream.ReadByte;
                s1 := Integer(ASourceStream.ReadByte) or ((b and $f0) shl 4);
                s2 := Integer(ASourceStream.ReadByte) or ((b and $0f) shl 8);

                Chunks[i, j * 2 + 0] := s1 - 2048;
                Chunks[i, j * 2 + 1] := s2 - 2048;
              end;
          else
            Assert(False, 'ChunkBitDepth not supported');
        end;

        FrameLength := ASourceStream.ReadDWord;

        SetLength(chunkIndex, ChannelCount);
        SetLength(chunkAttrs, ChannelCount);

        for i := 0 to FrameLength - 1 do
        begin
          for k := 0 to ChannelCount - 1 do
          begin
            w := ASourceStream.ReadWord;
            chunkIndex[k] := w and $fff;
            chunkAttrs[k] := w shr 12;
          end;

          for j := 0 to ChunkSize - 1 do
            for k := 0 to ChannelCount - 1 do
              memStream.WriteWord(Cardinal(CAttrLookup[chunkAttrs[k]] * Chunks[chunkIndex[k], j]) shr 15);
        end;
      end;

      memStream.Seek(0, soFromBeginning);
      ADestStream.Write(CreateWAVHeader(ChannelCount, 16, SampleRate, memStream.Size), SizeOf(TWavHeader));
      ADestStream.CopyFrom(memStream, memStream.Size);
    finally
      memStream.Free;
    end;
  end;


var
  gscFN, wavFN: String;
  inFS, outFS: TFileStream;
  gscFS: TMemoryStream;
begin
  if ParamCount = 0 then
  begin
    WriteLn('Usage: ', ExtractFileName(ParamStr(0)) + ' <source GSC file> [dest WAV file]');
    WriteLn;
    Exit;
  end;

  gscFN := ParamStr(1);
  if ParamCount > 1 then
    wavFN := ParamStr(2)
  else
    wavFN := ChangeFileExt(gscFN, '.wav');

  inFS := TFileStream.Create(gscFN, fmOpenRead or fmShareDenyNone);
  outFS := TFileStream.Create(wavFN, fmCreate or fmShareDenyWrite);
  gscFS := TMemoryStream.Create;
  try
    LZCompress(inFS, False, True, gscFS);
    gscFS.Position := 0;
    GSCUnpack(gscFS, outFS);
  finally
    inFS.Free;
    outFS.Free;
    gscFS.Free;
  end;
end.

