program decoder;

uses Types, SysUtils, Classes, Math, extern;

const
  CAttrMul = round((High(SmallInt) + 1) * (High(SmallInt) / 2047));
  CMaxAttenuation = 16;
  CVariableCodingHeaderSize = 2;
  CVariableCodingBlockSize = 3;

var
  CAttrLookup : array[Boolean {negative?}, 0 .. CMaxAttenuation - 1] of Integer;


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
    i, j, k, b, s1, s2, bitCount, variableCodingHeader: Integer;
    w: Word;
    chunkIndex: TIntegerDynArray;
    chunkNegative: TBooleanDynArray;
    StreamVersion, ChannelCount, ChunkBitDepth, ChunkSize, ChunkCount: Integer;
    FrameLength, SampleRate, ChunkBlend, AttenuationDivider: Integer;
    Chunks: TSmallIntDynArray2;
    Attenuations: TByteDynArray;
    memStream: TMemoryStream;
    law, lawAcc: Double;
    bits: Cardinal;

    function GetBits(ABitCount: Integer): Integer;
    begin
      Assert(ABitCount <= bitCount);
      Result := bits and ((1 shl ABitCount) - 1);
      bits := bits shr ABitCount;
      bitCount -= ABitCount;
    end;

  begin
    memStream := TMemoryStream.Create;
    try
      while ASourceStream.Position <> ASourceStream.Size do
      begin
        // parse header

        StreamVersion := ASourceStream.ReadByte;
        ChannelCount := ASourceStream.ReadByte;
        ChunkCount := ASourceStream.ReadWord and $1fff;
        ChunkBitDepth := ASourceStream.ReadByte;
        ChunkSize := ASourceStream.ReadByte;
        SampleRate := ASourceStream.ReadDWord;
        ChunkBlend := SampleRate shr 24;
        SampleRate := SampleRate and $ffffff;
        AttenuationDivider := ASourceStream.ReadWord;

        // compute attenuation law from AttenuationDivider

        law := 1.0 / AttenuationDivider;
        lawAcc := 1.0;
        for i := 0 to CMaxAttenuation - 1 do
        begin
          lawAcc += law * i;

          CAttrLookup[False, i] := round(CAttrMul / lawAcc);
          CAttrLookup[True, i] := -round(CAttrMul / lawAcc);
        end;

        if memStream.Position = 0 then
        begin
          writeln('StreamVersion = ', StreamVersion);
          writeln('SampleRate = ', SampleRate);
          writeln('ChannelCount = ', ChannelCount);
          writeln('ChunkBitDepth = ', ChunkBitDepth);
          writeln('ChunkSize = ', ChunkSize);
          writeln('ChunkBlend = ', ChunkBlend);
        end;

        Assert(ChunkBlend = 0, 'ChunkBlend not supported');

        SetLength(Chunks, ChunkCount, ChunkSize);
        SetLength(Attenuations, ChunkCount);

        // depack Attenuations

        for i := 0 to ChunkCount div 2 - 1 do
        begin
          b := ASourceStream.ReadByte;
          Attenuations[i * 2 + 0] := (b and $f0) shr 4;
          Attenuations[i * 2 + 1] := (b and $0f);
        end;

        if Odd(ChunkCount) then
        begin
          b := ASourceStream.ReadByte;
          Attenuations[ChunkCount - 1] := (b and $f0) shr 4;
        end;

        // depack Chunks

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

        // depack Frames

        FrameLength := ASourceStream.ReadDWord;

        SetLength(chunkIndex, ChannelCount);
        SetLength(chunkNegative, ChannelCount);

        bits := 0;
        bitCount := 0;
        variableCodingHeader := -1;
        for i := 0 to FrameLength - 1 do
        begin
          for k := 0 to ChannelCount - 1 do
          begin
            if (bitCount < 16) and (ASourceStream.Position < ASourceStream.Size) then
            begin
              w := ASourceStream.ReadWord;
              bits := bits or (w shl bitCount);
              bitCount += 16;
            end;

            chunkNegative[k] := GetBits(1) <> 0;

            if GetBits(1) <> 0 then // has new header?
              variableCodingHeader := GetBits(CVariableCodingHeaderSize);

            chunkIndex[k] := 0;
            for j := 0 to variableCodingHeader do
              chunkIndex[k] := (chunkIndex[k] shl CVariableCodingBlockSize) or GetBits(CVariableCodingBlockSize);
          end;

          for j := 0 to ChunkSize - 1 do
            for k := 0 to ChannelCount - 1 do
              memStream.WriteWord(Cardinal(CAttrLookup[chunkNegative[k], Attenuations[chunkIndex[k]]] * Chunks[chunkIndex[k], j]) shr 15);
        end;

        if bitCount >= 16 then // fixup potentially reading 1 spurious word
        begin
          ASourceStream.Seek(-2, soCurrent);
          bitCount -= 16;
        end;

        Assert(bitCount < 16);
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
  try
    GSCUnpack(inFS, outFS);
  finally
    outFS.Free;
    inFS.Free;
  end;
end.

