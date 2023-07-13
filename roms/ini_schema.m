function Nc = ini_schema(varargin)
%INI_SCHEMA Netcdf file schema for a ROMS initialization file
%
% This function creates a file schema for a ROMS initialization file.  The
% schema can be passed to ncwriteschema to create a file.  This is similar
% to the operations performed by c_initial.m in the myroms.org ROMS
% toolbox; however, returning the file schema instead of writing directly
% to file allows the user to make modifications as necessary (such as
% adding additional variables, adjusting attributes, etc.).
%
% Note that this function will only build the framework for the file; it 
% does not write to file itself, nor associate any specific data with the
% schema variables.
%
% Input variables (passed as parameter/value pairs):
%
%   spherical:  0 or 1, indicating Cartesian or spherical coordinates,
%               respectively [0]
%   
%   Vtransform: 1 or 2, indicating vertical transformation equation (1 =
%               original ROMS transformation, 2 = UCLA flavor
%               transformation) [1]
%
%   Lp:         number of *exterior* rho points in grid along the xi-axis
%               [1]
%
%   Mp:         number of *exterior* rho points in the grid along the
%               eta-axis [1]
%
%   N:          number of vertical levels [1]
%
%   NT:         number of active and passive travers [2]
%
%   NBL:        number of benthic layers [1]
%
%   format:     format for file, see nccreate for options ['classic']
%
%   reftime:    datetime, reference time used for time axis variable.
%               Units will be "seconds since reftime". [1900-01-01]
%
%   tname:      character array, name of time variable ['ocean_time']
%
%   vgroups:    string array or cell array of characters arrays, types of
%               variables to include in the file schema:
%               reference:  horizontal and vertical grid parameters
%               grid:       coordinate variables (lat/lon/x/y/s/h)
%               ocean:      ocean variables
%               ice:        ice variables
%               BEST_NPZ:   bio variables
%               BIO_COBALT: bio variables
%               BIO_BANAS:  bio variables
%
% Output variables:
%
%   Nc:         netCDF file schema.  See ncinfo, ncwriteschema for details.

% Copyright 2022 Kelly Kearney

p = inputParser;
p.addParameter('spherical', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '<=', 1}));
p.addParameter('Vtransform', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '<=', 2}));
p.addParameter('Lp', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('Mp', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('N',  1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('NT', 2, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('NBL', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('format', 'classic', @(x) validateattributes(x, {'string', 'char'}, {'scalartext'}));
p.addParameter('reftime', datetime(1900,1,1), @(x) validateattributes(x, {'datetime'}, {'scalar'}));
p.addParameter('tname', 'ocean_time', @(x) validateattributes(x, {'string', 'char'}, {'scalartext'}));
p.addParameter('vgroups', {'reference', 'grid', 'ocean'}, @(x) iscellstr(x) || isstring(x));

p.parse(varargin{:});
Opt = p.Results;

% Variable tables

V = romsvariablegroups(Opt);
 
% Parse variable groups

tf = ismember(Opt.vgroups, fieldnames(V));
if ~all(tf)
    str = join(fieldnames(V), ', ');
    error('Unrecognized variable group; options are %s', str);
end

V = cellfun(@(x) V.(x), Opt.vgroups, 'uni', 0);
V = cat(1, V{:});

% Initialize schema

Nc = ncschema_init(Opt.format);

% Global attributes

Nc = ncschema_addatts(Nc, 'type', 'INITIALIZATION file', ...
                          'history', sprintf('%s: %s', ...
                              datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                              'File schema created via ini_schema.m'));

% Add dimensions

Dims.r1 = {'s_rho'};
Dims.w1 = {'s_w'};
Dims.r2 = {'xi_rho' 'eta_rho', Opt.tname};
Dims.r3 = {'xi_rho' 'eta_rho' 's_rho', Opt.tname};
Dims.u2 = {'xi_u' 'eta_u', Opt.tname};
Dims.u3 = {'xi_u' 'eta_u' 's_rho', Opt.tname};
Dims.v2 = {'xi_v' 'eta_v', Opt.tname};
Dims.v3 = {'xi_v' 'eta_v' 's_rho', Opt.tname};
Dims.bn3 = {'xi_rho' 'eta_rho' 'benlayer', Opt.tname};
    
if Opt.spherical
    coords.r2 = 'lat_rho lon_rho';
    coords.u2 = 'lat_u lon_u';
    coords.v2 = 'lat_v lon_v';
    coords.r3 = 'lat_rho lon_rho s_rho';
    coords.u3 = 'lat_u lon_u s_rho';
    coords.v3 = 'lat_v lon_v s_rho';
    coords.bn3 = 'lat_rho lon_rho benlayer';
else
    coords.r2 = 'x_rho x_rho';
    coords.u2 = 'x_u x_u';
    coords.v2 = 'x_v x_v';
    coords.r3 = 'x_rho x_rho s_rho';
    coords.u3 = 'x_u x_u s_rho';
    coords.v3 = 'x_v x_v s_rho';
    coords.bn3 = 'x_rho y_rho benlayer';
end
coords = structfun(@(x) [x ' ' Opt.tname], coords, 'uni', 0);

Nc = ncschema_adddims(Nc, ...
    'xi_rho',       Opt.Lp,   false, ...
    'xi_u',         Opt.Lp-1, false, ...
    'xi_v',         Opt.Lp,   false, ...
    'xi_psi',       Opt.Lp-1, false, ...
    'eta_rho',      Opt.Mp,   false, ...
    'eta_u',        Opt.Mp,   false, ...
    'eta_v',        Opt.Mp-1, false, ...
    'eta_psi',      Opt.Mp-1, false, ...
    's_rho',        Opt.N,    false, ...
    's_w',          Opt.N+1,  false, ...
    'tracer',       Opt.NT,   false, ...
    'benlayer',     Opt.NBL,  false, ...
    Opt.tname,      1,        true);

dused = cellfun(@(x) Dims.(x), setdiff(unique(V.grid), 'nul'), 'uni', 0);
dused = unique(cat(2, dused{:}));

tf = ismember({Nc.Dimensions.Name}, [dused Opt.tname]);
Nc.Dimensions = Nc.Dimensions(tf);

% Add variables

Nc = ncschema_addvars(Nc, Opt.tname, ...
                      {Opt.tname}, ...
                      {'long_name', 'time since initialization', ...
                       'units', sprintf('seconds since %s', datestr(Opt.reftime, 'yyyy-mm-dd HH:MM:SS')), ...
                       'calendar', 'standard'}, ...
                      'double');


for iv = 1:height(V)
    
    % Build attribute array
    
    atts = {};
    if ~isempty(V.long{iv})
        atts = [atts 'long_name', V.long{iv}];
    end
    if ~isempty(V.units{iv})
        atts = [atts 'unit', V.units{iv}];
    end
    
    if isfield(V.grid{iv}, coords)
        atts = [atts 'coordinates', coord.(V.grid{iv})];
    end
    
    if ~ismember(V.grid{iv}, {'nul', 'r1', 'w1'})
        atts = [atts 'time', Opt.tname];
    end
    
    atts = [atts V.atts{iv}];
    
    % Add to schema
    
    if strcmp(V.grid{iv}, 'nul')
        dims = {};
    else
        dims = Dims.(V.grid{iv});
    end
    
    Nc = ncschema_addvars(Nc, V.short{iv}, dims, atts, V.class{iv});
    
end


