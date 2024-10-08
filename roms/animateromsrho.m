function h = animateromsrho(grdfile, bhis, this, varargin) %bsta, tsta, slat, slon)
%ANIMATEROMSRHO Plot interactive version of rho variables
%
% animateromsrho(grdfile, bhis, this)
% animateromsrho(Grd,     bhis, this)
% h = animateromsrho(Grd, bhis, this, p1, v1, ...)
%
% This function creates an interactive viewer that animates the spatial
% distriubtion of any ROMS rho variable.
%
% Input variables:
%
%   grdfile:    path to roms grid file
%
%   Grd:        ROMS grid file structure array (see ncstruct)
%
%   bhis:       nxi x neta x nt array of values to be animated
%
%   this:       nt x 1 array of datetimes corresponding to values to be
%               animated
%
% Optional input variables (passed as parameter/value pairs):
%
%   slat:       latitude of stations where timeseries should be
%               extracted/plotted
%
%   slon:       longitude of stations where timeseries should be
%               extracted/plotted
%
%   bsta:       nsta x ntsta array, timeseries corresponding to discrete
%               stations (if not included, values will be pulled from
%               nearest history grid cell)
%
%   tsta:       btsta x 1 array, datestimes corresponding to bsta
%               timeseries
%
%   

% Copyright 2016 Kelly Kearney

%--------------------
% Check data
%--------------------

% Parse inputs

p = inputParser;
p.addOptional('tsta', [], @(x) validateattributes(x, {'datetime'}, {'vector'}));
p.addOptional('bsta', [], @(x) validateattributes(x, {'numeric'}, {}));
p.addOptional('slat', [], @(x) validateattributes(x, {'numeric'}, {'<=', 90, '>=', -90}));
p.addOptional('slon', [], @(x) validateattributes(x, {'numeric'}, {'<=', 360, '>=', -180}));
p.parse(varargin{:});

tsta = p.Results.tsta;
bsta = p.Results.bsta;
slat = p.Results.slat;
slon = p.Results.slon;

% Grid file

if isstruct(grdfile)
    Grd = grdfile;
elseif ischar(grdfile) && exist(grdfile, 'file')
    Grd = ncstruct(grdfile);
else
    error('grdfile input must be ROMS gridfile structure or grid file name');
end

latlim = minmax(Grd.lat_rho);
lonlim = minmax(wrapTo360(Grd.lon_rho));

% Fill defaults where necessary

if isempty(slat)
    slat = mean(Grd.lat_rho(:));
end
if isempty(slon)
    slon = mean(Grd.lon_rho(:));
end

% Check data size

[nxi, neta] = size(Grd.lat_rho);
sz = size(bhis);
if ndims(bhis) > 3 || ~isequal(sz(1:2), [nxi neta])
    error('bhis must be nxi x neta x nt array');
end

% Extract "station" data if not provided separately

ell = referenceEllipsoid('earth');
distfun = @(ltln1, ltln2) distance(ltln1(:,1), ltln1(:,2), ltln2(:,1), ltln2(:,2), ell);
if isempty(bsta)
    idx = knnsearch([Grd.lat_rho(:) Grd.lon_rho(:)], [slat slon], 'Distance', distfun);
    [ixi,ieta] = ind2sub(size(Grd.lat_rho), idx);
    nstat = length(slat);
    bsta = nan(nstat, length(this));
    for ii = 1:nstat
        bsta(ii,:) = bhis(ixi(ii),ieta(ii),:);
    end
    tsta = this;
end

% Check for time duplicates

isrep = diff(this) == 0;
this = this(~isrep);
bhis = bhis(:,:,~isrep);

isrep = diff(tsta) == 0;
tsta = tsta(~isrep);
bsta = bsta(:,~isrep);

% Data limits

hislim = minmax(bhis(~isinf(bhis)));
stalim = minmax(bsta(~isinf(bsta)));
if stalim(1) == stalim(2)
    stalim = stalim(1) + [-1 1];
elseif any(isnan(stalim))
    stalim = [-1 1];
end

%--------------------
% Plot
%--------------------

h.fig = figure('color', 'w');

% Timeseries axis

h.tsax  = axes('position', [0.05 0.05 0.9 0.1], 'box', 'on');
hold on;
h.ln = plot(h.tsax, tsta, bsta);
h.ref = plot(tsta([1 1]), stalim, 'k');

% Map axis

h.mapax = axes('position', [0.05 0.15 0.9 0.8]);
worldmap(latlim, lonlim);

h.pc   = pcolorm(double(Grd.lat_psi), double(Grd.lon_psi), double(padend(bhis(2:end-1,2:end-1,1))));
% h.pc = pcolorm(Grd.lat_psi, Grd.lon_psi, padarray(bhis(2:end-1,2:end-1,1), [1 1], NaN, 'post'));
for ii = 1:length(slat)
    h.sta(ii) = plotm(slat(ii), slon(ii), 'o');
end
h.cb = colorbar('east');
setpos(h.cb, '# 0.6nz # 0.3nz');

set(h.mapax, 'clim', hislim);
lim = getm(h.mapax, 'maplatlimit');
setm(h.mapax, 'flinewidth', 1, 'mlabelparallel', lim(2));


set(h.sta, 'MarkerEdgeColor', 'w');
if isscalar(h.ln)
    set(h.sta, 'MarkerFaceColor', get(h.ln, 'color'));
else
    set(h.sta, {'MarkerFaceColor'}, get(h.ln, 'color'));
end

if verLessThan('matlab', '9.1') % R2016b
    set(h.tsax, 'ylim', stalim, 'xlim', datenum(minmax(tsta)));
else
    set(h.tsax, 'ylim', stalim, 'xlim', minmax(tsta));
    
end

set(h.tsax, 'ButtonDownFcn', @clickts);

% Buttons

pos = getpos(h.fig);
h.play = uicontrol('style', 'pushbutton', 'string', 'Play', ...
    'units', 'normalized', ...
    'position', [0.05 0.9 0.2 0.05]);
h.stop = uicontrol('style', 'pushbutton', 'string', 'Stop', ...
    'units', 'normalized', ...
    'position', [0.05 0.85 0.2 0.05]);
h.date = uicontrol('style', 'edit', 'string', char(this(1)), ...
    'units', 'normalized', ...
    'position', [0.05 0.8 0.2 0.05]);

t = timer('TimerFcn', @animatetime,...
        'ExecutionMode','fixedDelay',...
        'Period',0.1,...
        'TasksToExecute',Inf);
it = 1;

set(h.play, 'Callback', @(~,~) start(t));
set(h.stop, 'Callback', @(~,~) stop(t));

set(h.date, 'Callback', @(ht,~) updatetime(datetime(ht.String)));

set(h.fig, 'DeleteFcn', @closefig);


    function updatetime(t)
        
        [~,ihis] = min(abs(t - this));
        [~,ista] = min(abs(t - tsta));
        
        % set(h.pc, 'cdata', padarray(bhis(2:end-1,2:end-1,ihis), [1 1], NaN, 'post'));
        set(h.pc, 'cdata', double(padend(bhis(2:end-1,2:end-1,ihis))));
        if verLessThan('matlab', '9.1')
            set(h.ref, 'xdata', datenum(tsta([ista ista])));
        else
            set(h.ref, 'xdata', tsta([ista ista]));
        end
        if ~strcmp(h.date.String, char(t))
            h.date.String = char(t);
        end
        it = ihis;
        
    end

    function clickts(h, ~)
        xy = get(h, 'CurrentPoint');
        xlim = datenum(h.XLim);
        if xy(1,1) >= xlim(1) & xy(1,1) <= xlim(2) % is datenum
            updatetime(datetime(xy(1,1), 'convertfrom', 'datenum')); 
        else % is pixels?
            xnlim = h.XAxis.NumericLimits; % undocumented
            tnew = interp1(xnlim, h.XLim, xy(1,1));
            updatetime(tnew);
        end
    end

    function animatetime(src, ev)
        it = it + 1;
        tidx = mod(it, length(this));
        if tidx == 0
            tidx = length(this);
        end
        updatetime(this(tidx));
    end

    function closefig(src,ev)
        stop(t);
        delete(t);
    end

    function x = padend(x)
        x = [x nan(size(x,1),1); nan(1,size(x,2)+1)];
    end

end
