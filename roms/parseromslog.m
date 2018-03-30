function Log = parseromslog(file)
%PARSEROMSLOG Parse data from ROMS standard output
%
% The ROMS model prints several diagnostic variables (energy, basin volume,
% etc.) to standard output while it runs, interspersed with other
% information... setup details, input parameters, file reads and writes,
% etc.  This function parses some of that data.
%
% Input variables:
%
%   file:       name of text file where ROMS standard output was saved.
%
% Output variables:
%
%   Log:        structure with the following fields:
%
%               stime:      datatime, wall clock time for simulation start
%
%               etime:      datetime, wall clock time for simulation end
%
%               cputime:    table array displaying CPU time used by various
%                           ROMS calculations
%
%               timeref:    datetime, model reference time (i.e. date
%                           corresponding to STEP = 0)
%
%               dt:         duration, model time step
%
%               steps:      table array with time-stepping table data.
%                           Column names correspond to those in the file,
%                           with the exception of HH:MM:SS, which is
%                           changed to Time to  be consistent with Matlab
%                           variable name restrictions. A Date column is
%                           added that converts the Day+Time fields to a
%                           datetime based on the model reference time.

% Copyright 2017 Kelly Kearney

txt = fileread(file);
txt = regexp(txt, '\n', 'split')';

% Remove Node startup stuff

isnode = startsWith(txt, ' Node #');
txt = txt(~isnode);

% Parse clock start and end time

str = ' Model Input Parameters:  ROMS/TOMS version';
idx = find(startsWith(txt, str));

Log.stime = datetime(txt{idx+1}, 'inputformat', ...
    'eeee - MMMM d, yyyy - h:mm:ss a');

str = ' ROMS/TOMS: DONE...';
idx = find(startsWith(txt, str));

Log.etime = datetime(strtrim(strrep(txt{idx}, str, '')), ...
    'inputformat', 'eeee - MMMM d, yyyy - h:mm:ss a');

% Parse CPU time for various tasks

str1 = ' Nonlinear model elapsed time profile:';
str2 = ' Nonlinear model message Passage profile:';
str3 = ' All percentages are with respect to total time';

idx1 = find(startsWith(txt, str1));
idx2 = find(startsWith(txt, str2));
idx3 = find(startsWith(txt, str3));

txt1 = txt((idx1+2):(idx2-3));
txt2 = txt((idx2+2):(idx3-3));

tmp = regexp([txt1; txt2], '^([^\.]*)\s*\.*\s*([\d\.]*)', 'tokens', 'once');
tmp = cat(1, tmp{:});

% if ii == 1 && jj == 1
%     cpucats = tmp(:,1);
%     ncat = length(cpucats);
%     cputime = nan(nnode, nrst, ncat);
% end

Log.cputime = table(strtrim(tmp(:,1)), cellfun(@str2double, tmp(:,2)), 'VariableNames', {'Process', 'cpu_seconds'});

% A few time-related parameters

str = 'Reference time for units attribute (yyyymmdd.dd)';
idx = find(contains(txt, str));
tmp = regexp(strtrim(txt{idx}), '\s*', 'split');
time_ref = tmp{1};
yr = str2double(time_ref(1:4));
mn = str2double(time_ref(5:6));
dy = str2double(time_ref(7:end));
Log.timeref = datetime(yr, mn, dy);

str = 'Timestep size (s) for 3-D equations.';
idx = find(contains(txt, str));
tmp = regexp(strtrim(txt{idx}), '\s*', 'split');
Log.dt = seconds(str2double(tmp{1}));

% Main table

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
          
logtbl.Date = Log.timeref + days(logtbl.Day) + tdur;

Log.steps = logtbl;

