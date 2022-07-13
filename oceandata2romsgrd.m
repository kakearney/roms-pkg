function [vdata4, Rep] = oceandata2romsgrd(vdata, Grd, nlayer, varargin)
%DATA2ROMSGRD Interpolate/upsample 3D gridded data to ROMS grid
%
% out = oceandata2romsgrd(vdata, Grd, nlayer, lon, lat, depth)
% out = oceandata2romsgrd(vdata, Grd, nlayer, horiz, idx, ddist)
% [out, horiz, idx, ddist] = oceandata2romsgrd(...)
%
% This function performs interpolation/upsampling of a (usually)
% coarse-resolution, (usually) global gridded ocean data product to a ROMS
% grid, for use as initial conditions, boundary conditions, etc. 
%
% In most cases, the global data product will have NaNs/missing values near
% the coasts and near the bottom, where the coarse grid is land-masked but
% the finer-resolution ROMS grid requires data.  This function uses the
% following algorithm to fill these gaps and upsample the rest of the
% domain:
% 1) Fill bottom NaNs (below deepest depth in coarse grid) using nearest
%    neighbor from above.
% 2) Interpolate horizontally along each coarse-dataset depth surface, from
%    coarse grid to ROMS model grid, using lineaer interpolation (no
%    extrapolation), leaving NaNs in place.  This will fill interior points
%    but leave the NaN-blocks near coastlines.
% 3) Flood toward coastlines using inverse distance weighting with distance
%    accounting for barriers (i.e. the roms rho-grid land mask).  This
%    prevents interpolation through land.
% 4) Interpolate vertically, using linear interpolation and
%    nearest-neighbor extrapolation, from the fixed depth levels to
%    terrain-following ROMS depth levels. 
%
% Caveats:
%    - Interpolation steps treat lon/lat as cartesian coordinates for now.
%      Works fine if ROMS grid is approximate rectilinear, doesn't cross
%      any poles, and uses the same longitude conventions as the input
%      data. Should update to be more generic, though.
%    - The distance-with-barriers metric assumes the ROMS grid cells are
%      approximately square.  I think it should still work for more
%      irregular ROMS grids (given the inherent restrictions on the grids
%      for the hydrodynamic calculations), but I've never tested that. 
%
% Input variables:
%
%   vdata:  nlon x nlat x nz array, data to be interpolated 
%
%   Grd:    netcdf structure (see ncstruct) for ROMS grid file.
%
%   nlayer: number of depth levels for ROMS dataset
%
%   lon:    nlon x nlat array, or nlon x 1 vector, longitude coordinates of
%           gridded data 
%
%   lat:    nlon x nlat array, or nlat x 1 vector, latitude coordinates of
%           gridded data 
%
%   depth:  ndep x 1 array, depth coordinates of gridded data (m, positive
%           down) 
%
%   Rep:    See description below.  Can be passed as an input to reuse
%           previously-caclulated interpolation info with a new vdata array
%           (e.g., to save time when interpolating multiple datasets that
%           start on the same grid and target the same ROMS grid file)
%
% Output variables:
%
%   out:    nxi x neta x nlayer array, data interpolated to the ROMS grid
%
%   Rep:    structure of interpolation info associated with this specific
%           data grid + ROMS grid pairing, with the following fields:
% 
%           horiz:      griddedInterpolant object used to interpolate
%                       horizontally from the gridded dataset. 
%
%           ddist:      nfill x npt distance matrix of distances between
%                       fill-points (i.e. grid cells in the ROMS grid that
%                       require the coastal flooding algorithm) and
%                       gridded-data grid points
%
%           dataidx:    indices of the gridded-data horizontal grid points
%                       that are not land-masked at any depth 
%
%           needsfill:  logical mask indicating which grid cells in the
%                       ROMS grid require the coastal flooding
%
%           fvert:      griddedInterpolant object used to interpolate
%                       vertically from the gridded dataset

% Copyright 2022 Kelly Kearney

%------------------
% Setup
%------------------

% Check sizes/types

validateattributes(vdata, {'numeric'}, {'ndims',3});
[nlon,nlat,nz] = size(vdata);

validateattributes(Grd, {'struct'}, {});

validateattributes(nlayer, {'numeric'}, {'scalar', 'integer'});

% Check whether input includes lat/lon/z or previously-calculated
% interpolants

if isstruct(varargin{1}) && all(isfield(varargin{1}, {'horiz', 'ddist','dataidx','needsfill'}))
    reuse = true;
    
    horiz = varargin{1}.horiz;
    ddist = varargin{1}.ddist;
    dataidx = varargin{1}.dataidx;
    needsfill = varargin{1}.needsfill;
    fvert = varargin{1}.fvert;
    
else
    reuse = false;
    
    [lon,lat,depth] = deal(varargin{:});
    
    validateattributes(lon, {'double'}, {'>=',-180,'<=',360});
    validateattributes(lat, {'double'}, {'>=',-90,'<=',90});
    validateattributes(depth, {'numeric'}, {});
    
    if ~((isvector(lon) && isvector(lat)) || isequal(size(lon), size(lat)))
        error('lon and lat matrices should be vectors or ndgrid-ed matrices');
    end
    
    if isvector(lat) && isvector(lon)
        [lon, lat] = ndgrid(lon,lat);
    end
    
    if ~(isequal(size(lon), [nlon nlat]) && isequal(size(lat), [nlon nlat]))
        error('size mismatch between vdata and lon/lat');
    end
    
    if length(depth) ~= nz
        error('size mismatch between vdata and z');
    end
    
end
    
% ROMS grid size

[nxi, neta] = size(Grd.h);

%------------------
% Interpolate
%------------------

% Step 1: extend bottom data (in coarse-resolution) using nearest neighbor

fprintf('Extending vertically...\n');

bot = repmat(bottom(vdata), 1, 1, nz);
vdata(isnan(vdata)) = bot(isnan(vdata));
    
% Step 2: Upsample coarse-resolution layer to model grid horizontally

fprintf('Interpolating horizontally...\n');


vdata2 = nan(nxi, neta, nz);
if ~reuse
    horiz = griddedInterpolant({lon(:,1), lat(1,:)}, vdata(:,:,1));
end

for iz = 1:nz
    horiz.Values = vdata(:,:,iz);
    vdata2(:,:,iz) = horiz(Grd.lon_rho, Grd.lat_rho);
end

% Step 3a: Flood toward coastlines using inverse distance weighting

fprintf('Flooding coastlines...\n');

if ~reuse
    nlat = length(lat);
    nlon = length(lon);

    % Build interpolants to convert from lat/lon to xi/eta coordinates
    
    [xi,eta] = ndgrid(1:nxi, 1:neta);
    fxi = scatteredInterpolant(Grd.lat_rho(:), Grd.lon_rho(:), xi(:), 'nearest');
    feta = fxi; feta.Values = eta(:);

    % Determine which points in the coarse dataset are land-masked in the
    % coarse dataset but within the ROMS ocean grid

    [mln,mlt] = mask2poly(Grd.lon_psi, Grd.lat_psi, Grd.mask_rho(2:end-1,2:end-1));
    isin = inpolygon(lon, lat, mln, mlt);
    needsfill = isnan(vdata2(:,:,1)) & logical(Grd.mask_rho); 
    dataidx = find(~isnan(vdata(:,:,1)) & isin); % where there is data

    npt = length(dataidx);
    nfill = nnz(needsfill);
    ddist = nan(nfill, npt);

    % Convert coarse-grid data points to xi/eta space
    
    dxi = fxi(lat(dataidx), lon(dataidx));
    deta = feta(lat(dataidx), lon(dataidx));

    % Calculate distances from each data point to each ROMS grid point,
    % accounting for land barriers
    
    c = ConsoleProgressBar;
    c.setMaximum(npt);
    c.start();

    for ipt = 1:npt
        c.setValue(ipt);
        c.setText(sprintf('%d/%d\n', ipt, npt));
        t = traveltime(~Grd.mask_rho,  dxi(ipt), deta(ipt));
        ddist(:,ipt) = t(needsfill);
    end
    fprintf('\n');
end

% Step 3b: Apply inverse distance weighting

vdata3 = nan(nxi, neta, nz);
for iz = 1:nz
    tmp = vdata(:,:,iz);
    tmp = tmp(dataidx);
    u = idw(tmp', ddist, 5);

    tmp = vdata2(:,:,iz);
    tmp(needsfill) = u;

    vdata3(:,:,iz) = tmp;
end

% Step 4: Interpolate to depth surfaces

fprintf('Interpolating to depth levels...\n');

vdata4 = nan(nxi*neta,nlayer);
widx = find(Grd.mask_rho == 1);

if ~reuse
    fvert = griddedInterpolant(depth, ones(size(depth)), 'linear', 'nearest');
end

zr = calcromsz(Grd.h, 0, nlayer);

vtmp = cube2rect(vdata3, Grd.mask_rho==1);
ztmp = cube2rect(-zr, Grd.mask_rho==1);
for ig = 1:length(widx)
    fvert.Values = vtmp(:,ig);
    vdata4(widx(ig),:) = fvert(ztmp(:,ig));
end
vdata4 = reshape(vdata4, nxi, neta, nlayer);

% Step 5: fill closed-off areas with 0

for iz = 1:nlayer
    tmp = vdata4(:,:,iz);
    isn = Grd.mask_rho == 1 & isnan(tmp);
    tmp(isn) = 0;
    vdata4(:,:,iz) = tmp;
end

fprintf('Done\n');

% Save interpolation info for reuse

Rep.horiz = horiz;
Rep.ddist = ddist;
Rep.dataidx = dataidx;
Rep.needsfill = needsfill;
Rep.fvert = fvert;

% % Surface values spot check
% 
% h = plotgrid('size', [2 2]);
% axes(h.ax(1));
% worldmap(minmax(Grd.lat_rho), minmax(Grd.lon_rho));
% scatterm(lat(:), lon(:), [], reshape(vdata(:,:,1), [], 1), 'filled');
% bordersm('alaska','k'); bordersm('russia','k'); 
% 
% axes(h.ax(2));
% plotromsrho(Grd, vdata2(:,:,1));
% bordersm('alaska','k'); bordersm('russia','k'); 
% 
% axes(h.ax(3));
% plotromsrho(Grd, vdata3(:,:,1));
% bordersm('alaska','k'); bordersm('russia','k'); 
% 
% axes(h.ax(4));
% plotromsrho(Grd, vdata4(:,:,end));
% bordersm('alaska','k'); bordersm('russia','k'); 

%------------------------
% travel time subfunction
%------------------------

function H = traveltime(MM, ri, ci)
% This function was adapted from the following:
%
% Michael Kleder (2022). Shortest Path with Obstacle Avoidance (ver 1.3)
% (https://www.mathworks.com/matlabcentral/fileexchange/8625), MATLAB
% Central File Exchange.  
%
% It computes travel time (i.e. distance) between one grid point and all
% elements in a logical grid using finite element diffusion.  
%
% MM:   logical grid, 1 = barrier, 0 = clear
%
% ri:   row index of target point
%
% ci:   column index of target point

M=logical(MM);

speed = .5;
W=zeros(size(M));
W(ri,ci)=1;
H=zeros(size(M));
yespop = sum(W(:)>1);
nochangepop = 0;
itercount = 0;
while nochangepop < 20
    itercount = itercount + 1;
    diffconst = 0.1;
    W(1:end-1 ,:) = W(1:end-1 ,:) + diffconst * W(2:end,:);
    W(2:end,:) = W(2:end,:) + diffconst * W(1:end-1 ,:);
    W(:,1:end-1 ) = W(:,1:end-1 ) + diffconst * W(:,2:end );
    W(:,2:end) = W(:,2:end) + diffconst * W(:,1:end-1 );
    W(logical(M)) = 0;
    H(logical(W>1 & ~H))=itercount;
    yespopold = yespop;
    yespop = sum(W(:)>1);
    if abs(yespop-yespopold) < 1
        nochangepop = nochangepop + 1;
    end
end

H(H==0)=nan;


