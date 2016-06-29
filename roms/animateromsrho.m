function h = animateromsrho(grdfile, bhis, this, bsta, tsta, slat, slon)
%ANIMATEROMSRHO Plot interactive version of rho variables
%
% animateromsrho(grdfile, bhis, this, bsta, tsta, slat, slon)
% animateromsrho(Grd,     bhis, this, bsta, tsta, slat, slon)
%
% This function creates an interactive viewer that animates the spatial
% distriubtion of any ROMS rho variable.
%
% Input variables:
%


% Copyright 2016 Kelly Kearney

%--------------------
% Check data
%--------------------

% Grid file

if isstruct(grdfile)
    Grd = grdfile;
elseif ischar(grdfile) && exist(grdfile, 'file')
    Grd = ncreads(grdfile);
else
    error('grdfile input must be ROMS gridfile structure or grid file name');
end

latlim = minmax(Grd.lat_rho);
lonlim = minmax(wrapTo360(Grd.lon_rho));

% Check data size

[nxi, neta] = size(Grd.lat_rho);
sz = size(bhis);
if ndims(bhis) > 3 || ~isequal(sz(1:2), [nxi neta])
    error('bhis must be nxi x neta x nt array');
end

% Check for time duplicates

isrep = diff(this) == 0;
this = this(~isrep);
bhis = bhis(:,:,~isrep);

isrep = diff(tsta) == 0;
tsta = tsta(~isrep);
bsta = bsta(:,~isrep);

% Data limits

hislim = minmax(bhis);
stalim = minmax(bsta);
if stalim(1) == stalim(2)
    stalim = stalim(1) + [-1 1];
end

%--------------------
% Plot
%--------------------

h.fig = figure('color', 'w');

h.mapax = axes('position', [0.05 0.15 0.9 0.8]);
worldmap(latlim, lonlim);
h.pc = pcolorm(Grd.lat_psi, Grd.lon_psi, padarray(bhis(2:end-1,2:end-1,1), [1 1], NaN, 'post'));
for ii = 1:length(slat)
    h.sta(ii) = plotm(slat(ii), slon(ii), 'o');
end
h.cb = colorbar('east');
setpos(h.cb, '# 0.6nz # 0.3nz');

set(h.mapax, 'clim', hislim);
lim = getm(h.mapax, 'maplatlimit');
setm(h.mapax, 'flinewidth', 1, 'mlabelparallel', lim(2));

h.tsax  = axes('position', [0.05 0.05 0.9 0.1], 'box', 'on');
hold on;
h.ln = plot(h.tsax, tsta, bsta);
h.ref = plot(tsta([1 1]), stalim, 'k');

set(h.sta, 'MarkerEdgeColor', 'w');
set(h.sta, {'MarkerFaceColor'}, get(h.ln, 'color'));

set(h.tsax, 'ylim', stalim, 'xlim', datenum(minmax(tsta)));

set(h.tsax, 'ButtonDownFcn', @clickts);

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
        
        set(h.pc, 'cdata', padarray(bhis(2:end-1,2:end-1,ihis), [1 1], NaN, 'post'));
        set(h.ref, 'xdata', datenum(tsta([ista ista])));
        if ~strcmp(h.date.String, char(t))
            h.date.String = char(t);
        end
        it = ihis;
        
    end

    function clickts(h, ~)
        xy = get(h, 'CurrentPoint');
        updatetime(datetime(xy(1,1), 'convertfrom', 'datenum'));  
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

end
