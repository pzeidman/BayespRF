function data = Read_FreeSurfer(fname)
%data = Read_FreeSurfer(fname)
%
% Reads the contents of FreeSurfer ASCII files.
% The matrix is organized into five colums:
%   1:      vertex     vertex index number
%   2-4:    x,y,z      coordinates of vertices in world space
%   5:      values     contains overlay values etc
%
% If the file name has the extension '.label' the program reads a label
% (in which the first 2 lines are header information and thus ignored).
% If no file extension is given, the program assumes no headers.
%
% 02-04-2015:   Now permits loading csurf files with 6 columns (DSS)
%

fid = fopen(fname);
[p n e] = fileparts(fname);
try 
    if strcmpi(e,'.label')
        c = textscan(fid,repmat('%n',1,6),'headerlines',2);
    else
        c = textscan(fid,repmat('%n',1,5));
    end
catch
    data = NaN;
    return
end
fclose(fid);

data = cell2mat(c);

