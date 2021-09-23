function S = romsgeometryparams(Grd, varargin)
%ROMSGEOMETRYPARAMS Return parameters associated with ROMS grid geometry
%
% S = romsgeometryparams(Grd, varargin)
%
% This function provides a quick calculation of many of the geometry
% parameters related to a ROMS grid.
%
% Input variables:
%
%   Grd:            structure array, with fields corresponding to variables
%                   in a ROMS grid file (i.e. results of
%                   ncrstruct(grdname))  
%
% Optional input parameters:
%
%   spherical:      spherical grid switch, 1 for lat/lon spherical
%                   coordinates, 0 for x/y cartesian coordinates  
%                   [1]
%
%   N:              number of vertical levels at rho-points
%                   [10]
%
%   NT:             number of active and passive traver variables (usually
%                   2 for potential temperature and salinity)
%                   [2]
%
%   Vstretching:    stretching function to use
%                   1 = Song and Haidvogel (1994), original version
%                   2 = A. Shchepetkin (2005) UCLA-ROMS deprecated function
%                   3 = R. Geyer function for high bottom boundary layer
%                       resolution in relatively shallow applications
%                   4 = A. Shchepetkin (2010) UCLA-ROMS current function
%                   [1]
%
%   Vtransform:     transform function
%                   1 = original
%                   2 = UCLA-ROMS version
%                   [1]
%
%   theta_s:        S-coordinate surface control parameter
%                   [5]
%
%   theta_b:        S-coordinate bottom control parameter
%                   [0.4]
%
%   hc:             width of surface or bottom layer where higher vertical
%                   resolution is required (m) 
%                   [10]
%
%   tref:           datenumber, reference time for all time variables
%                   [datenum('1900/01/01')]
%                   
%
% Output variables:
%
%   S:              structure of geometry parameters, including all
%                   optional inputs plus:
%                   
%                   h:      Lr x Mr array, bathymetry (m, positive down) at
%                           rho points 
%
%                   [L/M][r/m/u/v]: scalars, number of rows (L) and columns
%                           (M) in the rho (r), interior-only rho (m), u
%                           (u) and v (v) grids of the domain  
%
%                   s_rho:  1 x N array, s-coordinate variable at rho
%                           points 
%
%                   Cs_r:   1 x N array, nondimensional, monotonic vertical
%                           stretching function, C(s), for rho points 
%
%                   s_w:    1 x N+1 array, s-coordinate variable at w
%                           points 
%
%                   Cs_w:   1 x N+1 array, nondimensional, monotonic
%                           vertical stretching function, C(w), for w
%                           points  
%
%                   z_r:    Lr x Mr x N array, depth (m, positive up) at
%                           rho points assuming zeta = 0
%
%                   z_u:    Lu x Mu x N array, depth (m, positive up) at
%                           u points assuming zeta = 0
%
%                   z_v:    Lv x Mv x N array, depth (m, positive up) at
%                           v points assuming zeta = 0
%
%                   z_w:    Lr x Mr x N+1 array, depth (m, positive up) at
%                           w points assuming zeta = 0
%
%                   Hz:     Lr x Mr x N array, thickness (m) of layers at
%                           rho points assuming zeta = 0

% Copyright 2017 Kelly Kearney

%-------------------------
% User inputs
%-------------------------

p = inputParser;
p.addParameter('spherical',   1,   @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative', '<=', 1}));
p.addParameter('N',          10,   @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}));
p.addParameter('NT',          2,   @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}));
p.addParameter('Vtransform',  1,   @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative', '<=', 2}));
p.addParameter('Vstretching', 1,   @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative', '<=', 4}));
p.addParameter('theta_s',     5,   @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('theta_b',     0.4, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('hc',         10,   @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('tref', datenum(1900,1,1), @(x) validateattributes(x, {'numeric'}, {'scalar'}));

p.parse(varargin{:});
S = p.Results;

%-------------------------
% Calculate depth-related 
% variables
%-------------------------

S.h = Grd.h;

[S.Lr, S.Mr] = size(Grd.h); % all rho points
S.Lm = S.Lr - 2;            % interior rho pointd
S.Mm = S.Mr - 2;
S.Lu = S.Lr - 1;            % u points 
S.Mu = S.Mr;    
S.Lv = S.Lr;                % v points 
S.Mv = S.Mr - 1;

%  Set vertical grid variables.

kgrid = 0;                                          % RHO-points

[S.s_rho, S.Cs_r] = stretching(S.Vstretching,                           ...
                             S.theta_s, S.theta_b, S.hc, S.N,         ...
                             kgrid, 0);

kgrid = 1;                                          % W-points

[S.s_w,   S.Cs_w] = stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N,         ...
                             kgrid, 0);

%  Compute ROMS model depths

ssh = zeros(size(Grd.h)); % Ignore free-surface contribution so 
                          % interpolation bounded below sea level

igrid = 1; % density points
S.z_r = set_depth(S.Vtransform, S.Vstretching,                        ...
                  S.theta_s, S.theta_b, S.hc, S.N,                    ...
                  igrid, S.h, ssh, 0);
      
igrid = 3; % u-velocity points
S.z_u = set_depth(S.Vtransform, S.Vstretching,                        ...
                  S.theta_s, S.theta_b, S.hc, S.N,                    ...
                  igrid, S.h, ssh, 0);

igrid = 4; % v-velocity points
S.z_v = set_depth(S.Vtransform, S.Vstretching,                        ...
                  S.theta_s, S.theta_b, S.hc, S.N,                    ...
                  igrid, S.h, ssh, 0);

%  Compute ROMS vertical level thicknesses (m).
      
igrid = 5; % w-velocity points
S.z_w = set_depth(S.Vtransform, S.Vstretching, ...
                  S.theta_s, S.theta_b, S.hc, S.N, ...
                  igrid, S.h, ssh, 0);

S.Hz = diff(S.z_w, 1, 3); 

