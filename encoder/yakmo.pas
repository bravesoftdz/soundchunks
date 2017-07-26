unit yakmo;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Types, Process, strutils;

procedure DoExternalKMeans(AFN: String; DesiredNbTiles, RestartCount: Integer; Normalize: Boolean; var XYC: TIntegerDynArray);

implementation

Const
  READ_BYTES = 65536; // not too small to avoid fragmentation when reading large files.

// helperfunction that does the bulk of the work.
// We need to also collect stderr output in order to avoid
// lock out if the stderr pipe is full.
function internalRuncommand(p:TProcess;var outputstring:string;
                            var stderrstring:string; var exitstatus:integer):integer;
var
    numbytes,bytesread,available : integer;
    outputlength, stderrlength : integer;
    stderrnumbytes,stderrbytesread : integer;
begin
  result:=-1;
  try
    try
    p.Options :=  [poUsePipes];
    bytesread:=0;
    outputlength:=0;
    stderrbytesread:=0;
    stderrlength:=0;
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
            Write(Copy(stderrstring, 1+StderrBytesRead, available));

            if StderrNumBytes > 0 then
              Inc(StderrBytesRead, StderrNumBytes);
          end
        else
          Sleep(100);
      end;
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
    internalRuncommand(Process, Output, ErrOut, i); // destroys Process

    DeleteFile(PChar(AFN));

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


end.

