function Log = parseromslog(file, varargin)
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
% Optional variables:
%
%   readsteps:  logical scalar, true to read time-stepping data (i.e. main
%               output table)
%
%   readcputime:logical scalar, true to read CPU useage table 
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
%
%               numcpu:     number of CPUs used (possibly mox-specific)


p = inputParser;
p.addParameter('readsteps', true);
p.addParameter('readcputime', true);
p.parse(varargin{:});

Opt = p.Results;

% Copyright 2017 Kelly Kearney

txt = fileread(file);
txt = regexp(txt, '\n', 'split')';

% Remove Node startup stuff

% isnode = startsWith(txt, ' Node #');
isnode = regexpfound(txt, '\s+Node\s+\#');

nodetxt = txt(isnode);
txt = txt(~isnode);

% Parse number of nodes

nodenum = regexp(nodetxt, 'Node\s*#([\s\d]*)', 'tokens', 'once');
nodenum = cellfun(@(x) str2num(x{1}), nodenum);
Log.numcpu = length(unique(nodenum));


% Parse clock start and end time

str = ' Model Input Parameters:  ROMS/TOMS version';
idx = find(startsWith(txt, str));

Log.stime = datetime(txt{idx+1}, 'inputformat', ...
    'eeee - MMMM d, yyyy - h:mm:ss a');

str = ' ROMS/TOMS: DONE...';
idx = find(startsWith(txt, str));

if ~isempty(idx)
    Log.etime = datetime(strtrim(strrep(txt{idx}, str, '')), ...
        'inputformat', 'eeee - MMMM d, yyyy - h:mm:ss a');
else
    Log.etime = NaT(0);
end

% Parse parameters

isparamtable = contains(txt, 'Parameters, Grid:');
idx = find(isparamtable);

isemp = @(x) isempty(x) || (ischar(x) & all(isspace(x)));

for ii = 1:length(idx)
    tmp = regexp(txt{idx(ii)}, '(\w*) Parameters, Grid: (\d*)', 'tokens', 'once');
    Log.paramtbl(ii).type = tmp{1};
    Log.paramtbl(ii).grid = str2double(tmp{2});
    
    idxstart = idx(ii)+3;
    idxend = find(cellfun(isemp, txt(idxstart:end)), 1)+idx(ii);
    tbltxt = txt(idxstart:idxend);
    
    tbltxt = strvcat(tbltxt{:});
    Log.paramtbl(ii).data = tbltxt;
    
    % phys grid 1: 11 18 54
    % bio  grid 1: 11 16 47
    % bio  grid 0: 11 16 71
    % ice  grid 1: 11 16 48 (extra line skip at top)
    % stat grid 1: 11 18 
    
    c1 = 11;
    switch Log.paramtbl(ii).type
        case {'Ice'}
            c2 = 16;
            tbltxt = tbltxt(2:end);  
        case {'Physical', 'Stations'}
            c2 = 18;
        case {'Biology'}
            c2 = 16;
        otherwise
            error('New type of parameter table');
    end
    val  = strtrim(cellstr(tbltxt(:,1:c1)));
    var  = strtrim(cellstr(tbltxt(:,c1+(1:c2))));
    desc = strtrim(cellstr(tbltxt(:,(c1+c2+1):size(tbltxt,2))));
   
    iscont = cellfun(@isempty,val) & cellfun(@isempty,var);
    idxc = find(iscont);
    for ic = 1:length(idxc)
        desc{idxc(ic)-1} = sprintf('%s %s', desc{idxc(ic)-1:idxc(ic)});
    end
    
    Log.paramtbl(ii).data = table(var(~iscont), val(~iscont), desc(~iscont), ...
        'variablenames', {'Var', 'Value', 'Description'});
end


% Parse CPU time for various tasks

if Opt.readcputime

    str1 = ' Nonlinear model elapsed time profile:';
    str2 = ' Nonlinear model message Passage profile:';
    str3 = ' All percentages are with respect to total time';

    idx1 = find(startsWith(txt, str1));
    idx2 = find(startsWith(txt, str2));
    idx3 = find(startsWith(txt, str3));

    if isempty(idx1) % Run terminated by outside process
        Log.cputime = [];
    else
        txt1 = txt((idx1+2):(idx2-3));
        txt2 = txt((idx2+2):(idx3-3));

        tmp = regexp([txt1; txt2], '^([^\.]*)\s*\.*\s*([\d\.]*)', 'tokens', 'once');
        tmp = cat(1, tmp{:});

        % if ii == 1 && jj == 1
        %     cpucats = tmp(:,1);
        %     ncat = length(cpucats);
        %     cputime = nan(nnode, nrst, ncat);
        % end

        Log.cputime = table(strtrim(tmp(:,1)), cellfun(@str2double, tmp(:,2)), ...
            'VariableNames', {'Process', 'cpu_seconds'});
    end
end

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

if Opt.readsteps

    idx = find(strncmp(txt, '   STEP   Day HH:MM:SS', 22));
    cols = txt{idx};

    txt = txt((idx+2):end);
    isemp = cellfun(@(x) isempty(strtrim(x)), txt);
    idx = find(isemp, 1);
    txt = txt(1:(idx-1));

    isdata = ~contains(txt, 'NaN') & regexpfound(txt, '^\s*\d');

    tmpfile = [tempname '.txt'];
    fid = fopen(tmpfile, 'wt');
    fprintf(fid, '%s\n', txt{isdata});
    fclose(fid);

    cols = regexp(strtrim(cols), '\s*', 'split');
    istime = strcmp(cols, 'HH:MM:SS');
    cols{istime} = 'Time';

    logtbl = readtable(tmpfile, 'ReadVariableNames', false);
    logtbl.Properties.VariableNames = cols;

    tdur = logtbl.Time;

    % tdur = duration(hour(datetime(logtbl.Time)), ...
    %               minute(datetime(logtbl.Time)), ...
    %               second(datetime(logtbl.Time)));

    logtbl.Date = Log.timeref + days(logtbl.Day) + tdur;

    Log.steps = logtbl;
end

