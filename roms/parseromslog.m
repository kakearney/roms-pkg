function logtbl = parseromslog(file, refdate)
%PARSEROMSLOG Parse the time-stepping table data from ROMS standard output
%
% The ROMS model prints several diagnostic variables (energy, basin volume,
% etc.) to standard output while it runs, interspersed with other
% information... setup details, input parameters, file reads and writes,
% etc.  This function parses only the time-stepping data.
%
% Input variables:
%
%   file:       name of text file where ROMS standard output was saved.
%
%   refdate:    date used as a reference time in the ROMS simulation, used
%               to convert time steps (in seconds) to dates.  If not
%               provided, the default assumes 1900-01-01.
%
% Output variables:
%
%   logtbl:     table array with time-stepping table data.  Column names
%               correspond to those in the file, with the exception of
%               HH:MM:SS, which is changed to Time to be consistent with
%               Matlab variable name restrictions.

% Copyright 2017 Kelly Kearney

if nargin < 2
    refdate = datetime(1900,1,1);
end

txt = fileread(file);
txt = regexp(txt, '\n', 'split');

idx = find(strncmp(txt, '   STEP   Day HH:MM:SS', 22));
cols = txt{idx};

txt = txt((idx+2):end);
isemp = cellfun(@(x) isempty(strtrim(x)), txt);
idx = find(isemp, 1);
txt = txt(1:(idx-1));

isdata = regexpfound(txt, '^\s*\d');

tmpfile = [tempname '.txt'];
fid = fopen(tmpfile, 'wt');
fprintf(fid, '%s\n', txt{isdata});
fclose(fid);

cols = regexp(strtrim(cols), '\s*', 'split');
istime = strcmp(cols, 'HH:MM:SS');
cols{istime} = 'Time';

logtbl = readtable(tmpfile, 'ReadVariableNames', false);
logtbl.Properties.VariableNames = cols;

tdur = duration(hour(datetime(logtbl.Time)), ...
              minute(datetime(logtbl.Time)), ...
              second(datetime(logtbl.Time)));
          
logtbl.Date = refdate + days(logtbl.Day) + tdur;

