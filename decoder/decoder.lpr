program decoder;

uses Types, SysUtils, Classes, Math, extern;

const
  CAttrMul = round((High(SmallInt) + 1) * (High(SmallInt) / 2047));
  CMaxAttenuation = 8;

var
  CAttrLookup : array[0 .. 15] of Integer = (
    round(CAttrMul / 1.0),
    round(CAttrMul / 1.2),
    round(CAttrMul / 1.6),
    round(CAttrMul / 2.2),
    round(CAttrMul / 3.0),
    round(CAttrMul / 4.0),
    round(CAttrMul / 5.2),
    round(CAttrMul / 6.6),
    -round(CAttrMul / 1.0),
    -round(CAttrMul / 1.2),
    -round(CAttrMul / 1.6),
    -round(CAttrMul / 2.2),
    -round(CAttrMul / 3.0),
    -round(CAttrMul / 4.0),
    -round(CAttrMul / 5.2),
    -round(CAttrMul / 6.6)
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
    StreamVersion, ChannelCount, ChunkBitDepth, ChunkSize, ChunkCount: Integer;
    FrameLength, SampleRate, ChunkBlend, AttenuationDivider: Integer;
    Chunks: TSmallIntDynArray2;
    memStream: TMemoryStream;
    law, lawAcc: Double;
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
        ChunkBlend := SampleRate shr 24;
        SampleRate := SampleRate and $ffffff;

        if StreamVersion > 1 then
        begin
          AttenuationDivider := ASourceStream.ReadWord;

          law := 1.0 / AttenuationDivider;
          lawAcc := 1.0;
          for i := 0 to CMaxAttenuation - 1 do
          begin
            lawAcc += law * i;

            CAttrLookup[i] := round(CAttrMul / lawAcc);
            CAttrLookup[i + CMaxAttenuation] := -round(CAttrMul / lawAcc);
          end;
        end;

        if memStream.Position = 0 then
        begin
          writeln('StreamVersion = ', StreamVersion);
          writeln('ChannelCount = ', ChannelCount);
          writeln('ChunkCount = ', ChunkCount);
          writeln('ChunkBitDepth = ', ChunkBitDepth);
          writeln('ChunkSize = ', ChunkSize);
          writeln('SampleRate = ', SampleRate);
          writeln('ChunkBlend = ', ChunkBlend);
        end;

        Assert(ChunkBlend = 0, 'ChunkBlend not supported');

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

