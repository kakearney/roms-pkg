function dt = ncdatelim(files, var)

if ischar(files)
    files = {files};
end
nfile = length(files);

t = nan(nfile,2);
for ii = 1:nfile
   
    I = ncinfo(files{ii}, var);
    if I.Size > 0
        t(ii,1) = ncread(files{ii}, var, 1, 1);
        t(ii,2) = ncread(files{ii}, var, I.Size, 1);
    end
end

tunit = ncreadatt(files{1}, var, 'units');
tparts = textscan(tunit, '%s since %D', 1);

switch lower(tparts{1}{1})
    case 'microseconds'
        dt = tparts{2} + seconds(t/1e6);
    case 'milliseconds'
        dt = tparts{2} + seconds(t/1000);
    case 'seconds'
        dt = tparts{2} + seconds(t);
    case 'minutes'
        dt = tparts{2} + minutes(t);
    case 'hours'
        dt = tparts{2} + hours(t);
    case {'days', 'day'}
        dt = tparts{2} + days(t);
    otherwise
        warning('Could not parse reference time');
        dt = [];
end
