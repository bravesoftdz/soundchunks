unit yakmo;

{$mode objfpc}{$H+}

interface

uses
  Windows, Classes, SysUtils, Types, Process, strutils, math;

type
  TDoubleDynArray2 = array of TDoubleDynArray;

procedure DoExternalKMeans(XY: TDoubleDynArray2;  DesiredNbTiles, RestartCount, Precision: Integer; AlternateMethod, PrintProgress: Boolean; var XYC: TIntegerDynArray);
function DoExternalEAQUAL(AFNRef, AFNTest: String; PrintStats, UseDIX: Boolean; BlockLength: Integer): Double;

implementation

Const
  READ_BYTES = 65536; // not too small to avoid fragmentation when reading large files.

// helperfunction that does the bulk of the work.
// We need to also collect stderr output in order to avoid
// lock out if the stderr pipe is full.
function internalRuncommand(p:TProcess;var outputstring:string;
                            var stderrstring:string; var exitstatus:integer; PrintOut: Boolean):integer;
var
    numbytes,bytesread,available : integer;
    outputlength, stderrlength : integer;
    stderrnumbytes,stderrbytesread, PrintLastPos, prp : integer;
begin
  result:=-1;
  try
    try
    p.Options :=  [poUsePipes];
    bytesread:=0;
    outputlength:=0;
    stderrbytesread:=0;
    stderrlength:=0;
    PrintLastPos:=1;
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

            // output to screen
            prp := Pos(#10, Copy(outputstring, PrintLastPos, bytesread - PrintLastPos + NumBytes));
            if PrintOut and (prp <> 0) then
            begin
              Write(Copy(outputstring, PrintLastPos, prp));
              PrintLastPos += prp;
            end;

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

            if StderrNumBytes > 0 then
              Inc(StderrBytesRead, StderrNumBytes);
          end
        else
          Sleep(20);
      end;

    if PrintOut then
      Write(Copy(stderrstring, PrintLastPos, StderrBytesRead - PrintLastPos));

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

procedure DoExternalKMeans(XY: TDoubleDynArray2; DesiredNbTiles, RestartCount, Precision: Integer; AlternateMethod,
  PrintProgress: Boolean; var XYC: TIntegerDynArray);
var
  i, j, Clu, Inp, st: Integer;
  InFN, Line, Output, ErrOut: String;
  SL, Shuffler: TStringList;
  Process: TProcess;
  OutputStream: TMemoryStream;
  Ratio: Double;
  pythonExe: array[0..MAX_PATH-1] of Char;
begin
  SL := TStringList.Create;
  Shuffler := TStringList.Create;
  OutputStream := TMemoryStream.Create;
  try
    for i := 0 to High(XY) do
    begin
      Line := IntToStr(i) + ' ';
      for j := 0 to High(XY[0]) do
        Line := Line + FloatToStr(XY[i, j]) + ' ';
      SL.Add(Line);
    end;

    InFN := GetTempFileName('', 'dataset-'+IntToStr(GetCurrentThreadId)+'.txt');
    SL.SaveToFile(InFN);
    SL.Clear;

    Process := TProcess.Create(nil);
    Process.CurrentDirectory := ExtractFilePath(ParamStr(0));

    if SearchPath(nil, 'python.exe', nil, MAX_PATH, pythonExe, nil) = 0 then
      pythonExe := 'python.exe';
    Process.Executable := pythonExe;

    for i := 0 to GetEnvironmentVariableCount - 1 do
      Process.Environment.Add(GetEnvironmentString(i));
    Process.Environment.Add('MKL_NUM_THREADS=1');
    Process.Environment.Add('NUMEXPR_NUM_THREADS=1');
    Process.Environment.Add('OMP_NUM_THREADS=1');

    Process.Parameters.Add('cluster.py');
    Process.Parameters.Add('-i "' + InFN + '" -n ' + IntToStr(DesiredNbTiles) + ' -t ' + FloatToStr(intpower(10.0, -Precision + 1)));
    if PrintProgress then
      Process.Parameters.Add('-d');
    if AlternateMethod then
      Process.Parameters.Add('-a');
    Process.ShowWindow := swoHIDE;
    Process.Priority := ppIdle;

    st := 0;
    internalRuncommand(Process, Output, ErrOut, st, PrintProgress); // destroys Process

    SL.LoadFromFile(InFN + '.membership');

    DeleteFile(PChar(InFN));
    DeleteFile(PChar(InFN + '.membership'));
    DeleteFile(PChar(InFN + '.cluster_centres'));

    SetLength(XYC, SL.Count);
    for i := 0 to SL.Count - 1 do
    begin
      Line := SL[i];
      XYC[i] := StrToIntDef(Line, -1);
    end;
  finally
    OutputStream.Free;
    Shuffler.Free;
    SL.Free;
  end;
end;


function DoExternalEAQUAL(AFNRef, AFNTest: String; PrintStats, UseDIX: Boolean; BlockLength: Integer): Double;
var
  i: Integer;
  Line, Output, ErrOut, SilFN: String;
  OutSL: TStringList;
  Process: TProcess;
  OutputStream: TMemoryStream;
begin
  SilFN := GetTempFileName('', 'silent-'+IntToStr(GetCurrentThreadId)+'.txt');

  Process := TProcess.Create(nil);
  OutSL := TStringList.Create;
  OutputStream := TMemoryStream.Create;
  try
    Process.CurrentDirectory := ExtractFilePath(ParamStr(0));
    Process.Executable := 'eaqual.exe';
    Process.Parameters.Add('-fref "' + AFNRef + '" -ftest "' + AFNTest + '"' + ifthen(BlockLength > 0, ' -blklen ' + IntToStr(BlockLength)));
    if not PrintStats then
      Process.Parameters.Add('-silent "' + SilFN + '"');
    Process.ShowWindow := swoHIDE;
    Process.Priority := ppIdle;

    i := 0;
    internalRuncommand(Process, Output, ErrOut, i, False); // destroys Process

    Result := -10.0;
    if PrintStats or not FileExists(SilFN) then
    begin
      OutSL.LineBreak := #13#10;
      OutSL.Text := Output;
      WriteLn(Output);
      WriteLn(ErrOut);

      for i := 0 to OutSL.Count - 1 do
      begin
        Line := OutSL[i];
        if (Pos('Resulting ODG:', Line) = 1) and not UseDIX or (Pos('Resulting DIX:', Line) = 1) and UseDIX then
        begin
          TryStrToFloat(RightStr(Line, Pos(#9, ReverseString(Line)) - 1), Result);
          Break;
        end;
      end;
    end
    else
    begin
      OutSL.LineBreak := #10;
      OutSL.LoadFromFile(SilFN);
      Line := OutSL[2];
      OutSL.Delimiter := #9;
      OutSL.DelimitedText := Line;
      TryStrToFloat(OutSL[Ord(UseDIX)], Result);
    end;

    DeleteFile(SilFN);

  finally
    OutputStream.Free;
    OutSL.Free;
  end;
end;

end.

