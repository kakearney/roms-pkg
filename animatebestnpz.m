function animatebestnpz(grdfile, hisfiles, stafile, var, iz, calc) %var, grdfile, hisfiles, stafile)
%ANIMATEBESTNPZ Plot interactive version of rho variables
%
% animatebestnpz(grdfile, hisfiles, stafile, var, iz, calc)
%
% This function creates an interactive viewer that animates the spatial
% distriubtion (from a slice of history file(s)) of a variable on a map
% while also plotting station data over time. 
%
% Input variables:
%
%   grdfile:    either the grid file or grid structure associated with the
%               ROMS run
%
%   hisfiles:   list of history files from the run
%
%   stafile:    station file from the run
%
%   var:        variable of interest
%
%   iz:         depth layer(s) of interest
%
%   calc:       function to apply to depth layers
%               'avg':  average values over depth
%               'sum':  integrate values over depth

% Copyright 2016 Kelly Kearney

% Get dimension info

Tmp = ncinfo(hisfiles{1});
isz = strcmp({Tmp.Dimensions.Name}, 'N');
nz = Tmp.Dimensions(isz).Length;

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

% Station and history dimensions

DimHis = collectromsdims(hisfiles, 'his', nz, Grd.h);
DimSta = collectromsdims(stafile, 'sta', nz);

DimHis.dz = diff(DimHis.zw, 1, 3);
DimSta.dz = diff(DimSta.zw, 1, 1);

% Read variable data from the indicated layers

isvar = strcmp({Tmp.Variables.Name}, var);
ndim = length(Tmp.Variables(isvar).Size);

tmp = cellfun(@(x) ncread(x, var), hisfiles, 'uni', 0);

if ndim == 3
    bhis = cat(3, tmp{:});       % nxi x neta x ntime 
    bsta = ncread(stafile, var); % nstation x ntime
elseif ndim == 4
    bhis = cat(4, tmp{:});       % nxi x neta x nz x ntime 
    bsta = ncread(stafile, var); % nz x nstation x ntime
    bsta(bsta > 1e35) = NaN;
    
    bhis = bhis(:,:,iz,:) .* DimHis.dz(:,:,iz,:);  
    bsta = bsta(iz,:,:) .* DimSta.dz(iz,:,:);
    
    switch calc
        case 'avg'
            bhis = sum(bhis,3)./sum(DimHis.dz(:,:,iz,:),3);
            bsta = sum(bsta,1)./sum(DimSta.dz(iz,:,:),1);   
        case 'sum'
            bhis = sum(bhis,3);
            bsta = sum(bsta,1);
        otherwise
            error('Unrecognized calc');
    end
    bhis = permute(bhis, [1 2 4 3]);
    bsta = permute(bsta, [2 3 1]);
      
end
this = datetime(1900,1,1) + DimHis.ocean_time/86400;
tsta = datetime(1900,1,1) + DimSta.ocean_time/86400;

hislim = minmax(bhis);
stalim = minmax(bsta);

% Set up plots

h.fig = figure('color', 'w');

h.mapax = axes('position', [0.05 0.15 0.9 0.8]);
worldmap(latlim, lonlim);
h.pc = pcolorm(Grd.lat_psi, Grd.lon_psi, padarray(bhis(2:end-1,2:end-1,1), [1 1], NaN, 'post'));
for ii = 1:length(DimSta.lat_rho)
    h.sta(ii) = plotm(DimSta.lat_rho(ii), DimSta.lon_rho(ii), 'o');
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





% [bhis, bhissum, bsta, bstasum] = deal(cell(ng,1));
% for ig = 1:ng
%     
%     tmp = cellfun(@(x) ncread(x, var), Files(ig).his, 'uni', 0);
%     bhis{ig} = cat(4, tmp{:});
%     
%     bhissum{ig} = sum(DimHis.dz .* bhis{ig}, 3);
%     
%     bsta{ig} = ncread(Files(ig).sta, var);
%     bstasum{ig} = sum(DimSta.dz .* bsta{ig}, 1);
%    
% end

% Plot

% cmap = cptcmap('Paired_12');
% 
% h = plotgrid('setup', cell(1,2), [], grps, 'sp', 0, 'mar', 0.05, 'mb', 0.1);
% for ii = 1:ng
%     axes(h.ax(ii));
%     worldmap(latlim, lonlim);
%     h.pc(ii) = pcolorm(Grd.lat_psi, Grd.lon_psi, padarray(bhissum{ig}(2:end-1,2:end-1,1,1), [1 1], NaN, 'post'));
%     
%     plotm(DimSta.lat_rho, DimSta.lon_rho, 'rx');
%     
% end
% 
% h.ax2 = axes('position', [0.05 0.05 0.9 0.1]);
% t = datetime(1900,1,1) + DimSta.ocean_time/86400;
% tmp = cat(1, bstasum{:});
% 
% h.ts = plot(t, reshape(tmp, [], size(tmp,3)));
% 
% set(h.ts, {'color'}, num2cell(cmap(1:6,:),2));
% 
% lim = minmax(tmp);
% set(h.ax, 'clim', lim);
% set(h.ax2, 'ylim', lim, 'xlim', minmax(datenum(t)));
% 
% hold(h.ax2, 'on');
% h.ref = plot([t(1) t(1)], lim, 'color', 'k', 'parent', h.ax2);
% 
% 
% 
% 
% 
% for it = 1:size(bhissum{1},4)
%     for ig = 1:2
%         set(h.pc(ig), 'cdata', padarray(bhissum{ig}(2:end-1,2:end-1,1,it), [1 1], NaN, 'post'));
%     end
%     
%     [~, imin] = min(abs(DimHis.ocean_time(it) - DimSta.ocean_time));
%     set(h.ref, 'xdata', datenum(t([imin imin])));
% 
%     pause;
% end
