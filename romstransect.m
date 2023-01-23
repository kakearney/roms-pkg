function h = romstransect(Grd, val, varargin)
%ROMSTRANSECT Extract a transect from a ROMS rho-grid variable
%
% h = romstransect(Grd, val)
% h = romstransect(Grd, val, p1, v1, ...)
%
% Input variables:
%
%   Grd:    ROMS grid structure
%
%   val:    ROMS rho-grid variable (nxi x neta x srho)
%
% Optional input variables:
%
%   xi:     xi indices of transects
%
%   eta:    eta indices of transects
%
%   zeta:   zeta value, scalar or nxi x neta grid
%
%   expand: true to expand values to equal sizes for pcolor plotting



[nxi, neta] = size(Grd.h);

p = inputParser;
p.addParameter('xi', 0, @(x) validateattributes(x, {'numeric'}, {'integer', '<=', nxi}));
p.addParameter('eta', 0, @(x) validateattributes(x, {'numeric'}, {'integer', '<=', neta}));
p.addParameter('zeta', 0, @(x) isscalar(x) || isequal(size(x), [nxi neta]));
p.addParameter('expand', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.parse(varargin{:});

Opt = p.Results;

nlayer = size(val,3);

Opt.xi = Opt.xi(Opt.xi>0);
Opt.eta = Opt.eta(Opt.eta>0);


for ixi = 1:length(Opt.xi)
        
    v = permute(val(Opt.xi(ixi),:,:), [3 2 1]);
    lt = Grd.lat_rho(Opt.xi(ixi),:);
    ln = Grd.lon_rho(Opt.xi(ixi),:);

    if isscalar(Opt.zeta)
        zeta = Opt.zeta;
    else
        zeta = Opt.zeta(Opt.xi(ixi),:);
    end

    [~,zw] = calcromsz(Grd.h(Opt.xi(ixi),:), zeta, nlayer);
    zw = permute(zw, [3 2 1]);

    x = distance(lt(1), ln(1), lt, ln, referenceEllipsoid('earth', 'km'));

    if Opt.expand
        x = repmat(x, nlayer+1, 1);
        lt = repmat(lt, nlayer+1, 1);
        ln = repmat(ln, nlayer+1, 1);
        v = cat(1, v, nan(1, size(x,2)));
    end
    
    h.xi(ixi).idx = Opt.xi(ixi);
    h.xi(ixi).x = x;
    h.xi(ixi).lt = lt;
    h.xi(ixi).ln = ln;
    h.xi(ixi).z = zw;
    h.xi(ixi).v = v;

end

for ieta = 1:length(Opt.eta)
        
    v = permute(val(:,Opt.eta(ieta),:), [3 1 2]);
    lt = Grd.lat_rho(:,Opt.eta(ieta))';
    ln = Grd.lon_rho(:,Opt.eta(ieta))';

    if isscalar(Opt.zeta)
        zeta = Opt.zeta;
    else
        zeta = Opt.zeta(:,Opt.eta(ieta));
    end

    [~,zw] = calcromsz(Grd.h(:,Opt.eta(ieta)), zeta, nlayer);
    zw = permute(zw, [3 1 2]);

    x = distance(lt(1), ln(1), lt, ln, referenceEllipsoid('earth', 'km'));

    if Opt.expand
        x = repmat(x, nlayer+1, 1);
        lt = repmat(lt, nlayer+1, 1);
        ln = repmat(ln, nlayer+1, 1);
        v = cat(1, v, nan(1, size(x,2)));
    end
    
    h.eta(ieta).idx = Opt.eta(ieta);
    h.eta(ieta).x = x;
    h.eta(ieta).lt = lt;
    h.eta(ieta).ln = ln;
    h.eta(ieta).z = zw;
    h.eta(ieta).v = v;

end