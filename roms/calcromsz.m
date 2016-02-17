function [zr, zw] = calcromsz(h, zeta, nlayer, varargin)
%CALCROMSZ Calculate ROMS depth (z) values at rho and omega points
%
% [zr, zw] = calcromsz(h, zeta, nlayer, p1, v1)
%
% Input variables:
%
%   h:              array of bottom depth values (m, positive down).
%                   Typically a 2D xi  x eta grid, but can be any
%                   dimensions
%
%   zeta:           array of free surface height, (m, negative down).
%                   Dimensions must be [size(h) nt] (nt can be 1, i.e. same
%                   size as h).
%
%   nlayer:         number of depth layers in grid
%
% Optional input variables (passed as parameter/value pairs or structure)
%
%   Vstretching:    stretching function to use [1]
%                   1 = Song and Haidvogel (1994), original version
%                   2 = A. Shchepetkin (2005) UCLA-ROMS deprecated function
%                   3 = R. Geyer function for high bottom boundary layer
%                       resolution in relatively shallow applications
%                   4 = A. Shchepetkin (2010) UCLA-ROMS current function
%
%   Vtransform:     transform function [1]
%                   1 = original
%                   2 = UCLA-ROMS version
%
%   theta_s:        S-coordinate surface control parameter [5]
%
%   theta_b:        S-coordinate bottom control parameter [0.4]
%
%   hc:             width of surface or bottom layer where higher vertical
%                   resolution is required (m) [10]
%
% Output variables:
%
%   zr:             array size [size(h) nlayer nt], depth at rho points,
%                   i.e. center of each layer (m) 
%
%   zw:             array size [size(h) nlayer+1 nt], depth at omega
%                   points, i.e. edges of each layer (m)   

% Copyright 2015 Kelly Kearney

p = inputParser;
p.addParameter('Vtransform',  1,   @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '<=', 2}));
p.addParameter('Vstretching', 1,   @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '<=', 4}));
p.addParameter('theta_s',     5,   @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('theta_b',     0.4, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('hc',         10,   @(x) validateattributes(x, {'numeric'}, {'scalar'}));

p.KeepUnmatched = true;

p.parse(varargin{:});
Opt = p.Results;

% TODO: add check to make sure zeta is same size as h or one dimension
% longer

% Calculate s and C

[sr, Cr] = stretching(Opt.Vstretching, Opt.theta_s, Opt.theta_b, Opt.hc, ...
                      nlayer, 0, false);
[sw, Cw] = stretching(Opt.Vstretching, Opt.theta_s, Opt.theta_b, Opt.hc, ...
                      nlayer, 1, false);

% Calculate z

hsz = size(h);
zsz = size(zeta);

if isequal(hsz, zsz)
    nt = 1;
else
    nt = zsz(end);
end

npt = numel(h);
h = h(:);

hr = h * ones(1,nlayer);
hw = h * ones(1, nlayer+1);

zeta = permute(reshape(zeta, [], nt), [1 3 2]);
zetar = repmat(zeta, [1 nlayer 1]);
zetaw = repmat(zeta, [1 nlayer+1 1]);

sw = ones(npt,1) * sw;
sr = ones(npt,1) * sr;
Cw = ones(npt,1) * Cw;
Cr = ones(npt,1) * Cr;

switch Opt.Vtransform
    case 1
        z0r = Opt.hc .* (sr - Cr) + hr .* Cr;
        z0w = Opt.hc .* (sw - Cw) + hw .* Cw;
        
        if nt > 1
            z0r = repmat(z0r, [1 1 nt]);
            z0w = repmat(z0w, [1 1 nt]);
        end
        
        zr = z0r + zetar .* (1 + z0r./repmat(hr, [1 1 nt]));
        zw = z0w + zetaw .* (1 + z0w./repmat(hw, [1 1 nt]));
        
    case 2
        z0r = (Opt.hc .* sr + Cr .* hr)./(hr + Opt.hc);
        z0w = (Opt.hc .* sw + Cw .* hw)./(hw + Opt.hc);
        
        if nt > 1
            z0r = repmat(z0r, [1 1 nt]);
            z0w = repmat(z0w, [1 1 nt]);
        end
        
        zr = zetar + (zetar + repmat(hr, [1 1 nt])) .* z0r;
        zw = zetaw + (zetaw + repmat(hw, [1 1 nt])) .* z0w;
end

zr = reshape(zr, [hsz nlayer nt]);
zw = reshape(zw, [hsz nlayer+1 nt]);
