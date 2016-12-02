function [t, tdate] = readoceantime(files)
%READOCEANTIME Read ROMS ocean_time variable
%
% [t, tdate] = readoceantime(files)
%
% This function reads the ocean_time variable form a ROMS output file, and
% also converts the values to datetimes, based on the ocean_time units
% attribute.
%
% Input variables:
%
%   files:  cell array of files
%
% Output variables:
%
%   t:      n x 1 array, ocean_time values
%
%   tdate:  n x 1 array, ocean_time values converted to datetimes

% Copyright 2016 Kelly Kearney

t = cellfun(@(x) ncread(x, 'ocean_time'), files, 'uni', 0);
t = cat(1, t{:});

tunit = ncreadatt(files{1}, 'ocean_time', 'units');
tparts = textscan(tunit, '%s since %D', 1);
switch lower(tparts{1}{1})
    case 'seconds'
        tdate = tparts{2} + seconds(t);
    case 'hours'
        tdate = tparts{2} + hours(t);
    case 'days'
        tdate = tparts{2} + days(t);
    otherwise
        warning('Could not parse reference time');
        tdate = [];
end