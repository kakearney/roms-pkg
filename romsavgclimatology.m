function [Out, Aligned] = romsavgclimatology(yavg, tavg, varargin)
%ROMSAVGCLIMATOLOGY Calculate climatology from ROMS averages output
%
% Clim = romsavgclimatology(yavg, tavg)
% Clim = romsavgclimatology(yavg, tavg, p1, v1, ...)
% [Clim, Aligned] = romsavgclimatology(yavg, tavg, p1, v1, ...)
%
% This function calculates a climatology from an input timeseries, with the
% assumption that each time step represents a binned average (centered on
% the indicated time value) rather than a single snapshot.
%
% Input variables:
%
%   yavg:       n x 1 vector, values to be climatologically averaged
%
%   tavg:       n x 1 vector of datetimes, mid-point times corresponding to
%               values
%
% Optional input values:
%
%   nbin:       Number of time bins to use in output climatology (365~daily,
%               52~weekly).
%               [52]
%
%   prctile:    vector of length np, percentile values to return, 0-100.
%               [25 50 75]
%
%   realign:    logical scalar, if true, also returns the timeseries
%               realigned into a nyear x nbin array corresponding to the
%               bins used to calculate the climatology.
%
% Output variables:
%
%   Clim:        structure with the following fields:
%
%               doy:    1 x nbin, day of year corresponding to midpoint of
%                       each climatological bin, *relative to pivot date*
%               mean:   1 x nbin, climatological mean
%               std:    1 x nbin, climatological standard deviation
%               prc:    np x nbin, climatological percentiles
%
%   Aligned:    structure with the following fields:
%
%               y:      nyr x nbin array, input y-values reweighted from
%                       their original time bins to those used for the
%                       climatology 
%           
%               t:      nyr x nbin+1 array, datetime values corresponding
%                       to the edges of each climotological time bin        

% Copyright 2023 Kelly Kearney

p = inputParser;
p.addParameter('nbin', 52, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
p.addParameter('prctile', [25 50 75], @(x) validateattributes(x, {'numeric'}, {'vector', '>=', 0, '<=', 100}));
p.addParameter('realign', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('pivotmonth', 1, @(x) validateattributes(x, {'numeric'}, {'vector', 'integer', '<=', 12, '>=', 1}));
p.addParameter('pivotday', 1, @(x) validateattributes(x, {'numeric'}, {'vector', 'integer', '<=', 31, '>=', 1}));
p.parse(varargin{:});

Opt = p.Results;

% Remove NaNs

isn = isnan(yavg);
yavg = yavg(~isn);
tavg = tavg(~isn);

dtavg = mean(diff(tavg));

% Determine years in timeseries (adjusting as necessary to account for a
% shifted pivot date)

tlim = minmax(tavg);
if tlim(1) < datetime(year(tlim(1)),Opt.pivotmonth,Opt.pivotday)
    t1 = datetime(year(tlim(1))-1,Opt.pivotmonth,Opt.pivotday);
else
    t1 = datetime(year(tlim(1)),Opt.pivotmonth,Opt.pivotday);
end

if tlim(2) <= datetime(year(tlim(2)),Opt.pivotmonth,Opt.pivotday)
    t2 = datetime(year(tlim(2)),Opt.pivotmonth,Opt.pivotday);
else
    t2 = datetime(year(tlim(2))+1,Opt.pivotmonth,Opt.pivotday);
end

% yrs = unique(year(tavg));
% 
% if min(tavg) < datetime(yrs(1), Opt.pivotmonth, Opt.pivotday)
%     yrs = [yrs(1)-1; yrs];
% end
% if max(tavg) > datetime(yrs(end), Opt.pivotmonth, Opt.pivotday)
%     yrs = [yrs; yrs(end)+1];
% end

yrs = year(t1):(year(t2)-1);

nyr = length(yrs);
tbin = NaT(length(yrs), Opt.nbin+1);

% Define bins evenly-spaced intervals across each year.  Leap years
% will have slightly larger intervals.

for ii = 1:length(yrs)
    tbin(ii,:) = linspace(datetime(yrs(ii),  Opt.pivotmonth,Opt.pivotday), ...
                          datetime(yrs(ii)+1,Opt.pivotmonth,Opt.pivotday), Opt.nbin+1);
end

% Consolidate list of bin edges, removing the duplicate pivot dates

tbinv = [reshape(tbin(1:end-1,1:end-1)', 1, []), tbin(end,:)];

% Calculate how much of each average-time-step interval falls into each
% climatology bin.

tbinv = datenum(tbinv); % b/c Matlab won't implictly expand datetimes
thi = datenum(tavg + dtavg/2);
tlo = datenum(tavg - dtavg/2);

dtperbin = max(min(thi, tbinv(2:end)) - max(tlo, tbinv(1:end-1)), 0);

% If requested, calculate realigned data

if Opt.realign
    weightperbin = dtperbin./sum(dtperbin,1);
    Aligned.y = reshape(sum(yavg.*weightperbin, 1), Opt.nbin, length(yrs))';
    Aligned.t = tbin;
end

% Aggregate single-week bins across years.  We now have a matrix showing
% how much each data point contributes to each climatology bin.

weeknum = repmat(1:Opt.nbin, 1, nyr);
dtperbin = splitapply(@sum, dtperbin', weeknum')';

%---------------
% Weighted stats
%---------------

weight = dtperbin./sum(dtperbin,1);
% weight(isnan(weight)) = 0; % TODO: Need to consider how to handle NaNs

% Time

Out.doy = mean(days((tbin(:,1:end-1) + diff(tbin,1,2)/2) - tbin(:,1)),1);

% Mean 

Out.mean = sum(yavg .* weight, 1);

% Standard deviation

Out.std = nan(1, Opt.nbin);
for ii = 1:Opt.nbin
    if ~all(isnan(weight(:,ii)))
        Out.std(ii) = std(yavg, weight(:,ii));
    end
end

% Percentiles

nprc = length(Opt.prctile);

Out.prc = nan(nprc, Opt.nbin);
for ii = 1:Opt.nbin
    if ~all(isnan(weight(:,ii)))
        Out.prc(:,ii) = wprctile(yavg, Opt.prctile, weight(:,ii), 5);
    end
end





