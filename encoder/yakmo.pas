unit yakmo;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Types, Process, strutils;

procedure DoExternalKMeans(AFN: String; DesiredNbTiles, RestartCount: Integer; Normalize: Boolean; var XYC: TIntegerDynArray);
function DoExternalEAQUAL(AFNRef, AFNTest: String; PrintStats, UseDIX: Boolean; BlockLength: Integer): Double;
function DoExternalEAQUALMulti(AFNRef, AFNTest: String; UseDIX: Boolean; BlockCount, BlockLength: Integer): TDoubleDynArray;

implementation

Const
  READ_BYTES = 65536; // not too small to avoid fragmentation when reading large files.

// helperfunction that does the bulk of the work.
// We need to also collect stderr output in order to avoid
// lock out if the stderr pipe is full.
function internalRuncommand(p:TProcess;var outputstring:string;
                            var stderrstring:string; var exitstatus:integer; PrintErr: Boolean):integer;
var
    numbytes,bytesread,available : integer;
    outputlength, stderrlength : integer;
    stderrnumbytes,stderrbytesread, PrintErrLastPos, perp : integer;
begin
  result:=-1;
  try
    try
    p.Options :=  [poUsePipes];
    bytesread:=0;
    outputlength:=0;
    stderrbytesread:=0;
    stderrlength:=0;
    PrintErrLastPos:=1;
    p.Execute;
    while p.Running do
      begin
        // Only call ReadFromStream if Data from corresponding stream
        // is already available, otherwise, on  linux, the read call
        // is blocking, and thus it is not possible to be sure to handle
        // big data amounts bboth on output and stderr pipes. PM.
        available:=P.Output.NumBytesAvailable;
        if  available > 0 then
          begin
            if (BytesRead + available > outputlength) then
              begin
                outputlength:=BytesRead + READ_BYTES;
                Setlength(outputstring,outputlength);
              end;
            NumBytes := p.Output.Read(outputstring[1+bytesread], available);
            if NumBytes > 0 then
              Inc(BytesRead, NumBytes);
          end
        // The check for assigned(P.stderr) is mainly here so that
        // if we use poStderrToOutput in p.Options, we do not access invalid memory.
        else if assigned(P.stderr) and (P.StdErr.NumBytesAvailable > 0) then
          begin
            available:=P.StdErr.NumBytesAvailable;
            if (StderrBytesRead + available > stderrlength) then
              begin
                stderrlength:=StderrBytesRead + READ_BYTES;
                Setlength(stderrstring,stderrlength);
              end;
            StderrNumBytes := p.StdErr.Read(stderrstring[1+StderrBytesRead], available);

            // output stderr to screen
            perp := Pos(#10, Copy(stderrstring, PrintErrLastPos, StderrBytesRead - PrintErrLastPos + stderrnumbytes));
            if PrintErr and (perp <> 0) then
            begin
              Write(Copy(stderrstring, PrintErrLastPos, perp));
              PrintErrLastPos += perp;
            end;

            if StderrNumBytes > 0 then
              Inc(StderrBytesRead, StderrNumBytes);
          end
        else
          Sleep(40);
      end;

    if PrintErr then
      Write(Copy(stderrstring, PrintErrLastPos, StderrBytesRead - PrintErrLastPos));

    // Get left output after end of execution
    available:=P.Output.NumBytesAvailable;
    while available > 0 do
      begin
        if (BytesRead + available > outputlength) then
          begin
            outputlength:=BytesRead + READ_BYTES;
            Setlength(outputstring,outputlength);
          end;
        NumBytes := p.Output.Read(outputstring[1+bytesread], available);
        if NumBytes > 0 then
          Inc(BytesRead, NumBytes);
        available:=P.Output.NumBytesAvailable;
      end;
    setlength(outputstring,BytesRead);
    while assigned(P.stderr) and (P.Stderr.NumBytesAvailable > 0) do
      begin
        available:=P.Stderr.NumBytesAvailable;
        if (StderrBytesRead + available > stderrlength) then
          begin
            stderrlength:=StderrBytesRead + READ_BYTES;
            Setlength(stderrstring,stderrlength);
          end;
        StderrNumBytes := p.StdErr.Read(stderrstring[1+StderrBytesRead], available);
        if StderrNumBytes > 0 then
          Inc(StderrBytesRead, StderrNumBytes);
      end;
    setlength(stderrstring,StderrBytesRead);
    exitstatus:=p.exitstatus;
    result:=0; // we came to here, document that.
    except
      on e : Exception do
         begin
           result:=1;
           setlength(outputstring,BytesRead);
         end;
     end;
  finally
    p.free;
  end;
end;

procedure DoExternalKMeans(AFN: String; DesiredNbTiles, RestartCount: Integer; Normalize: Boolean; var XYC: TIntegerDynArray);
var
  i, Clu, Inp: Integer;
  Line, Output, ErrOut: String;
  OutSL: TStringList;
  Process: TProcess;
  OutputStream: TMemoryStream;
begin
  Process := TProcess.Create(nil);
  OutSL := TStringList.Create;
  OutputStream := TMemoryStream.Create;
  try
    Process.CurrentDirectory := ExtractFilePath(ParamStr(0));
    Process.Executable := 'yakmo.exe';
    Process.Parameters.Add('"' + AFN + '" - - -O 2 -k ' + IntToStr(DesiredNbTiles) + ' -m ' + IntToStr(RestartCount) + IfThen(Normalize, ' -n'));
    Process.ShowWindow := swoHIDE;
    Process.Priority := ppIdle;

    i := 0;
    internalRuncommand(Process, Output, ErrOut, i, True); // destroys Process

    DeleteFile(PChar(AFN));

    if (Pos(#10, Output) <> Pos(#13#10, Output) + 1) then
      OutSL.LineBreak := #10;

    OutSL.Text := Output;

    for i := 0 to OutSL.Count - 1 do
    begin
      Line := OutSL[i];

      Inp := StrToInt(Copy(Line, 1, Pos(' ', Line) - 1));
      Clu := StrToIntDef(RightStr(Line, Pos(' ', ReverseString(Line)) - 1), 0);

      XYC[Inp] := Clu;
    end;
  finally
    OutputStream.Free;
    OutSL.Free;
  end;
end;


function DoExternalEAQUAL(AFNRef, AFNTest: String; PrintStats, UseDIX: Boolean; BlockLength: Integer): Double;
var
  i, Clu, Inp: Integer;
  Line, Output, ErrOut: String;
  OutSL: TStringList;
  Process: TProcess;
  OutputStream: TMemoryStream;
begin
  Process := TProcess.Create(nil);
  OutSL := TStringList.Create;
  OutputStream := TMemoryStream.Create;
  try
    Process.CurrentDirectory := ExtractFilePath(ParamStr(0));
    Process.Executable := 'eaqual.exe';
    Process.Parameters.Add('-fref "' + AFNRef + '" -ftest "' + AFNTest + '" -blklen ' + IntToStr(BlockLength));
    Process.ShowWindow := swoHIDE;
    Process.Priority := ppIdle;

    i := 0;
    internalRuncommand(Process, Output, ErrOut, i, False); // destroys Process

    OutSL.LineBreak := #13#10;
    OutSL.Text := Output;

    if PrintStats then
      WriteLn(Output);

    for i := 0 to OutSL.Count - 1 do
    begin
      Line := OutSL[i];
      if (Pos('Resulting ODG:', Line) = 1) and not UseDIX or (Pos('Resulting DIX:', Line) = 1) and UseDIX then
      begin
        Result := -StrToFloatDef(RightStr(Line, Pos(#9, ReverseString(Line)) - 1), 0);
        Break;
      end;
    end;
  finally
    OutputStream.Free;
    OutSL.Free;
  end;
end;

function DoExternalEAQUALMulti(AFNRef, AFNTest: String; UseDIX: Boolean; BlockCount, BlockLength: Integer): TDoubleDynArray;
var
  i: Integer;
  Line, Output, ErrOut, OutFN: String;
  OutSL: TStringList;
  Process: TProcess;
  OutputStream: TMemoryStream;
begin
  Process := TProcess.Create(nil);
  OutSL := TStringList.Create;
  OutputStream := TMemoryStream.Create;
  try
    OutFN := GetTempFileName('', 'out-'+IntToStr(GetCurrentThreadId)+'.txt');

    Process.CurrentDirectory := ExtractFilePath(ParamStr(0));
    Process.Executable := 'eaqual.exe';
    Process.Parameters.Add('-fref "' + AFNRef + '" -ftest "' + AFNTest + '" -blockout ' + IfThen(UseDIX, 'DI', 'ODG') + ' "' + OutFN + '" -blklen ' + IntToStr(BlockLength));
    Process.ShowWindow := swoHIDE;
    Process.Priority := ppIdle;

    Assert(FileExists(AFNRef));
    Assert(FileExists(AFNTest));
    Assert(not FileExists(OutFN));

    i := 0;
    internalRuncommand(Process, Output, ErrOut, i, False); // destroys Process

    Assert(FileExists(OutFN));

    OutSL.LineBreak := #13#10;
    OutSL.LoadFromFile(OutFN);

    SetLength(Result, BlockCount);
    for i := 0 to BlockCount - 1 do
    begin
      Line := OutSL[i * 2 + 2];
      Result[i] := -StrToFloatDef(Line, 0.0);
    end;

    DeleteFile(OutFN);

  finally
    OutputStream.Free;
    OutSL.Free;
  end;
end;

end.

