function [vdata4, Rep, I] = oceandata2romsgrd(vdata, Grd, nlayer, varargin)
%OCEANDATA2ROMSGRD Interpolate/upsample 3D gridded data to ROMS grid
%
% [out, Rep] = oceandata2romsgrd(vdata, Grd, nlayer, p1, v1)
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
% 1) Interpolate horizontally along each coarse-dataset depth surface, from
%    coarse grid to ROMS model grid, using lineaer interpolation (no
%    extrapolation), leaving NaNs in place.  This will fill interior points
%    but leave the NaN-blocks near coastlines.
% 2) Flood toward coastlines using inverse distance weighting with distance
%    accounting for barriers (i.e. the roms rho-grid land mask).  This
%    prevents interpolation through land.
% 3) Flood vertically via nearest neighbor to fill any missing data at
%    depth from above
% 4) Interpolate vertically, using linear interpolation and
%    nearest-neighbor extrapolation, from the fixed depth levels to
%    terrain-following ROMS depth levels. 
%
% Caveats:
%    - Interpolation steps treat lon/lat as cartesian coordinates for now.
%      Works fine if ROMS grid is approximately rectilinear, doesn't cross
%      any poles, and uses the same longitude conventions as the input
%      data. Should update to be more generic, though.
%    - The distance-with-barriers metric assumes the ROMS grid cells are
%      approximately square.  I think it should still work for more
%      irregular ROMS grids (given the inherent restrictions on the grids
%      for the hydrodynamic calculations), but I've never tested that. 
%
%
% Required input variables:
%
%   vdata:  nlon x nlat x nz array, data to be interpolated 
%
%   Grd:    netcdf structure (see ncstruct) for ROMS grid file.
%
%   nlayer: number of depth levels for ROMS dataset
%
% Input variables (semi-optional, passed as parameter/value pairs):
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
%   Rep:    See description in output.  Can be passed as an input to reuse
%           previously-caclulated interpolation info with a new vdata
%           array.  This saves significant time when interpolating multiple
%           datasets that start on the same grid and target the same ROMS
%           grid.
%
%   The following combinations of the above options are allowed:
%
%   2D interpolation: 
%       oceandata2romsgrd(..., 'lat', lt, 'lon', ln);  
%   3D interpolation: 
%       oceandata2romsgrd(..., 'lat', lt, 'lon', ln, 'depth', z); 
%   3D interpolation reusing a 2D interpolant:
%       oceandata2romsgrd(..., 'Rep', Rep, 'depth', z); 
%   2/3D interpolation reusing a interpolant of same:
%       oceandata2romsgrd(..., 'Rep', Rep); 
%
% Optional input variables (passed as parameter/value pairs):
%
%   grid:   1 x 2 string array or cell array of strings, indicating the
%           horizontal and vertical ROMS grid, respectively, for output 
%           dataset. 
%           horizontal options: 'rho', 'u', 'v', 'psi'
%           vertical options: 'rho', 'w'
%           The following combinations correspond to ROMS grid types:
%           ["rho", "rho"] = rho-grid (default)
%           ["u",   "rho"] = u-grid
%           ["v",   "rho"] = v-grid
%           ["psi", "rho"] = psi-grid
%           ["rho", "w"  ] = w-grid
%
%   zeta:   xi_rho x eta_rho, free surface height (m) corresponding to the
%           ROMS rho-grid (values will be interpolated to other grids as
%           necessary).  If not included, assumed to be 0 everywhere.
%
%   nedge:  size of square filter used to identify near-coastline grid
%           cells in the coarse-resolution dataset.  Only these will be
%           used to build the interpolants. [3]
%
%   verbose:logical scalar, true to print progress statements [true]
%
%   gpoly:  polyshape object for grid water mask.  If not included, this
%           will be calculated using mask2poly.m.  However, this function
%           may fail on grids if water-edge vertices connect to more than
%           two line segments; in these instances, it is safer to
%           pre-calculate the water masking polygon to check for a
%           well-formed polygon.
%
%   fillspeed:  diffusion constant (0 < fillspeed < 1) used for coastal
%           flooding algorithm. [0.1]
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
%           flood:      nz x 1 structure related to horizontal flooding of 
%                       coastal grid cells in each depth level
%
%                       ddist:      nfill x npt distance matrix of
%                                   distances between  fill-points (i.e.
%                                   grid cells in the ROMS grid that
%                                   require the coastal flooding algorithm)
%                                   and gridded-data grid points   
%
%                       dataidx:    indices of the gridded-data horizontal 
%                                   grid points that are not land-masked at
%                                   any depth
%
%                       needsfill:  logical mask indicating which grid
%                                   cells in the ROMS grid require the
%                                   coastal flooding 
%
%           fvert:      griddedInterpolant object used to interpolate
%                       vertically from the gridded dataset

% Copyright 2022 Kelly Kearney

%------------------
% Setup
%------------------

p = inputParser;

p.addParameter('grid', ["rho" "rho"], @(x) validateattributes(x, {'string', 'cell'}, {}));
p.addParameter('lon', [], @(x) validateattributes(x, {'double'}, {'>=',-180,'<=',360}));
p.addParameter('lat', [], @(x) validateattributes(x, {'double'}, {'>=',-90, '<=',90}));
p.addParameter('depth', [], @(x) validateattributes(x, {'double'}, {}));
p.addParameter('Rep', struct, @(x) validateattributes(x, {'struct'}, {}));
p.addParameter('zeta', 0, @(x) validateattributes(x, {'double'}, {}));
p.addParameter('nedge', 3, @(x) validateattributes(x, {'numeric'}, {'integer', 'odd'}));
p.addParameter('verbose', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('gpoly', [], @(x) validateattributes(x, {'polyshape'}, {}));
p.addParameter('fillspeed', 0.1, @(x) validateattributes(x, {'numeric'}, {'>', 0, '<', 1}));
p.addParameter('vinterp', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.parse(varargin{:});
Opt = p.Results;

% Check sizes/types

validateattributes(Grd, {'struct'}, {});

validateattributes(nlayer, {'numeric'}, {'scalar', 'integer'});
is2D = nlayer == 1;

if is2D
    validateattributes(vdata, {'numeric'}, {'ndims',2});
else
    validateattributes(vdata, {'numeric'}, {'ndims',3});
end
[nlon,nlat,nz] = size(vdata);

grid(1) = validatestring(Opt.grid(1), ["rho", "psi", "u", "v"], 'oceandata2romsgrd', 'grid(1)', 3);
grid(2) = validatestring(Opt.grid(2), ["rho", "w"],             'oceandata2romsgrd', 'grid(2)', 3);

% Check whether input includes lat/lon/z or previously-calculated
% interpolants

if isstruct(Opt.Rep) && all(isfield(Opt.Rep, {'horiz', 'flood'}))
    reuse = true;
    
    horiz = Opt.Rep.horiz;
    flood = Opt.Rep.flood;
%     ddist = Opt.Rep.ddist;
%     dataidx = Opt.Rep.dataidx;
%     needsfill = Opt.Rep.needsfill;

    if ~is2D
        if isfield(Opt.Rep, 'fvert')
            fvert = Opt.Rep.fvert;
        else
            fvert = [];
            depth = Opt.depth;
            validateattributes(depth, {'numeric'}, {});
            if length(depth) ~= nz
                error('size mismatch between vdata and z');
            end
        end
    end
else
    reuse = false;
    
    lon = Opt.lon;
    lat = Opt.lat;
    depth = Opt.depth;
    
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
    
    if ~is2D
        if length(depth) ~= nz
            error('size mismatch between vdata and z');
        end
    end
    
end

% Output grid

Out.lat = Grd.("lat_" + grid(1));
Out.lon = Grd.("lon_" + grid(1));
Out.mask = Grd.("mask_" + grid(1));

if ~is2D
    [zr, zw] = calcromsz(Grd.h, Opt.zeta, nlayer);
    switch grid(2)
        case "rho"
            ztmp = zr;
        case "w"
            ztmp = zw;
    end
    switch grid(1)
        case "rho"
            Out.z = ztmp;
        case "psi"
            Out.z = (ztmp(2:end,1:end-1,:) + ztmp(2:end,2:end,:) + ztmp(1:end-1,1:end-1,:) + ztmp(1:end-1,2:end,:))./4;
        case "u"
            Out.z = (ztmp(1:end-1,:,:) + ztmp(2:end,:,:))./2;
        case "v"
            Out.z = (ztmp(:,1:end-1,:) + ztmp(:,2:end,:))./2;
    end
end

% ROMS grid size

[nxi, neta] = size(Out.mask);

%------------------
% Interpolate
%------------------

% Step 1: extend bottom data (in coarse-resolution) using nearest neighbor

% if ~is2D
%     fprintf('Extending vertically...\n');
%     
% %     % Check if there are any nearby grid cells with data at this depth that
% %     % might be able to flood these.
% %    
% %     for iz = 1:nz
% %         tmp = movmean(vdata(:,:,iz), 3, 1, 'omitnan', 'Endpoints', 'shrink');
% %         tmp = movmean(tmp, 3, 2, 'omitnan', 'Endpoints', 'shrink');
% %     end
% 
%     bot = repmat(bottom(vdata), 1, 1, nz);
%     vdata(isnan(vdata)) = bot(isnan(vdata));
% end
    
% Step 2: Upsample coarse-resolution layer to model grid horizontally

if Opt.verbose
    fprintf('Interpolating horizontally...\n');
end
    
vdata2 = nan(nxi, neta, nz);
if ~reuse
    horiz = griddedInterpolant({lon(:,1), lat(1,:)}, vdata(:,:,1));
end

for iz = 1:nz
    horiz.Values = vdata(:,:,iz);
    vdata2(:,:,iz) = horiz(Out.lon, Out.lat);
end

% Step 3a: Flood toward coastlines using inverse distance weighting

if Opt.verbose
    fprintf('Flooding coastlines...\n');
end

if ~reuse
    nlat = length(lat);
    nlon = length(lon);

    % Build interpolants to convert from lat/lon to xi/eta coordinates
    
    [xi,eta] = ndgrid(1:nxi, 1:neta);
    fxi = scatteredInterpolant(Out.lat(:), Out.lon(:), xi(:), 'nearest');
    feta = fxi; feta.Values = eta(:);

    % Determine which points in the coarse dataset are land-masked in the
    % coarse dataset but within the ROMS ocean grid
    
    if isempty(Opt.gpoly)
        [mln,mlt] = mask2poly(Grd.lon_psi, Grd.lat_psi, Grd.mask_rho(2:end-1,2:end-1));
    else
        [mln, mlt] = boundary(Opt.gpoly);
    end
    dataingrid = inpolygon(lon, lat, mln, mlt);
    
    se = strel('square', Opt.nedge);
    for iz = 1:nz
        hasdata = ~isnan(vdata(:,:,iz));
        if is2D
            iswater = Out.mask == 1;
        else
            iswater = -Out.z(:,:,1) > depth(iz) & Out.mask == 1;
        end
        flood(iz).needsfill = isnan(vdata2(:,:,iz)) & iswater;
        
        isdataedge = hasdata & ~imerode(hasdata,se); % limit points used to those near the data edge
        flood(iz).dataidx = find(isdataedge & dataingrid);
       
        npt = length(flood(iz).dataidx);
        nfill = nnz(flood(iz).needsfill);
        flood(iz).ddist = nan(nfill, npt);
        
        dxi  = fxi( lat(flood(iz).dataidx), lon(flood(iz).dataidx));
        deta = feta(lat(flood(iz).dataidx), lon(flood(iz).dataidx));
        
%         test = traveltimemulti(~iswater, dxi, deta);
        if Opt.verbose
            c = ConsoleProgressBar;
            c.setMaximum(npt);
            c.start();
        end
%         test2 = nan(size(test));
        for ipt = 1:npt
            if Opt.verbose
                c.setValue(ipt);
                c.setText(sprintf('Depth %d/%d, %d/%d\n', iz, nz, ipt, npt));
            end
            t = traveltime(~iswater,  dxi(ipt), deta(ipt), Opt.fillspeed);
%             test2(:,:,ipt) = t;
            flood(iz).ddist(:,ipt) = t(flood(iz).needsfill);
        end
        
        if Opt.verbose
            fprintf('\n');
        end
    end
end

% Step 3b: Apply inverse distance weighting

vdata3 = nan(nxi, neta, nz);
for iz = 1:nz
    tmp = vdata(:,:,iz);
    tmp = tmp(flood(iz).dataidx);
    u = idw(tmp', flood(iz).ddist, 5);

    tmp = vdata2(:,:,iz);
    tmp(flood(iz).needsfill) = u;

    vdata3(:,:,iz) = tmp;
end

% Step 4a: Fill vertically

if ~is2D

    bot = repmat(bottom(vdata3), 1, 1, nz);
    vdata3(isnan(vdata3)) = bot(isnan(vdata3));
    
end


% Step 4: Interpolate to depth surfaces

if is2D || ~Opt.vinterp
    vdata4 = vdata3;
else
    
    if Opt.verbose
        fprintf('Interpolating to depth levels...\n');
    end

    vdata4 = nan(nlayer,nxi*neta);
    widx = Out.mask == 1;
% 
%     if ~reuse || isempty(fvert)
%         fvert = griddedInterpolant(depth, ones(size(depth)), 'linear', 'nearest');
%     end
% 
%     vtmp = cube2rect(vdata3, Out.mask==1);
%     ztmp = cube2rect(-Out.z, Out.mask==1);
%     for ig = 1:length(widx)
%         fvert.Values = vtmp(:,ig);
%         vdata4(widx(ig),:) = fvert(ztmp(:,ig));
%     end
%     vdata4 = reshape(vdata4, nxi, neta, nlayer);
    
    vtmp = cube2rect(vdata3, Out.mask==1);
    ztmp = cube2rect(-Out.z, Out.mask==1);

    z1 = repmat(fvert.GridVectors{1}, 1, size(vtmp,2));
    ztmp = min(max(ztmp, repmat(z1(1,:),30,1)), repmat(z1(end,:),30,1));

    vdata4(:,widx) = fastinterpcol(z1, vtmp, ztmp(end:-1:1,:));
    vdata4 = vdata4(end:-1:1,:)';
    vdata4 = reshape(vdata4, nxi, neta, nlayer);
    
    % Step 5: fill closed-off areas with 0

    for iz = 1:nlayer
        tmp = vdata4(:,:,iz);
        isn = Out.mask == 1 & isnan(tmp);
        tmp(isn) = 0;
        vdata4(:,:,iz) = tmp;
    end

    if Opt.verbose
        fprintf('Done\n');
    end
end

% Save interpolation info for reuse

Rep.horiz = horiz;
Rep.flood = flood;
% Rep.ddist = ddist;
% Rep.dataidx = dataidx;
% Rep.needsfill = needsfill;
if ~is2D
    Rep.fvert = fvert;
end

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

function H = traveltime(MM, ri, ci, diffconst)
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

W=zeros(size(M));
W(ri,ci)=1;
H=zeros(size(M));
yespop = sum(W(:)>1);
nochangepop = 0;
itercount = 0;
while nochangepop < 20
    itercount = itercount + 1;
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

function H = traveltimemulti(MM, ri, ci)

M=logical(MM);
ni = length(ri);
M3 = repmat(M, 1, 1, ni);

speed = .5;
W=zeros([size(M), ni]);
for ii = 1:ni
    W(ri(ii),ci(ii),ii)=1;
end
H=zeros([size(M), ni]);

counter = zeros(1,ni);

% yespop = sum(reshape(W, [],ni) > 1, 1);

nochangepop = 0;
itercount = 0;

underthresh = true(1,ni);

diffconst = 0.1;

while any(counter < 20)
    mask = counter < 20;
    
    itercount = itercount + 1;
%     fprintf('%d\n', itercount);

    check1 = W > 1;
    
    % Diffuse
    
    W(1:end-1,:,mask) = W(1:end-1 ,:,mask) + diffconst * W(2:end,:,mask);
    W(2:end,  :,mask) = W(2:end,:,mask) + diffconst * W(1:end-1 ,:,mask);
    W(:,1:end-1,mask) = W(:,1:end-1,mask) + diffconst * W(:,2:end,mask);
    W(:,2:end,  mask) = W(:,2:end,mask) + diffconst * W(:,1:end-1,mask);
    W(M3) = 0;
    
    % If value over 1 in cell and wasn't previously, mark with iteration
    
    H(W>1 & ~H & permute(mask,[1 3 2]))=itercount;
    
    % Did at least one grid cell newly rise over threshold?
    
    check2 = W > 1;
    
    isnew = check2 & ~check1;
    
    counter = counter + ~any(reshape(isnew, [], ni), 1);
    
%     yespopold = yespop;
%     yespop = sum(reshape(W, [],ni) > 1, 1);
%     
%     underthresh = (yespop - yespopold) < 1;
%     if any(underthresh)
%         nochangepop = nochangepop + 1;
%     end
%     if abs(yespop-yespopold) < 1
%         nochangepop = nochangepop + 1;
%     en

end
H(H==0) = NaN;




