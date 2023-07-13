function Nc = bry_schema(varargin)
%BRY_SCHEMA Build ROMS boundary file schema
%
% This function creates a file schema for a ROMS boundary file.  The
% schema can be passed to ncwriteschema to create a file.  This is similar
% to the operations performed by c_boundary.m in the myroms.org ROMS
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
%   boundary:   string array or cell array of character arrays, boundary
%               directions to include.  Must be a subset of the default
%               [["west" "east" "south" "north"]] 
%
%   tname:      character array, name of time variable ['bry_time']
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

% Copyright 2023 Kelly Kearney

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
p.addParameter('boundary', ["west" "east" "south" "north"], @(x) validateattributes(x, {'string','cell'}, {}));
p.addParameter('tname', 'bry_time', @(x) validateattributes(x, {'string', 'char'}, {'scalartext'}));
p.addParameter('vgroups', {'reference', 'grid', 'ocean'}, @(x) iscellstr(x) || isstring(x));

p.parse(varargin{:});
Opt = p.Results;

if iscell(Opt.boundary)
    if ~iscellstr(Opt.boundary)
        error('boundary must be either string or cell array of character arrays');
    end
    Opt.boundary = string(Opt.boundary);
end
for ib = 1:length(Opt.boundary)
    Opt.boundary(ib) = validatestring(Opt.boundary(ib), ["west" "east" "south" "north"], 'bry_schema', 'boundary');
end
Opt.boundary = reshape(Opt.boundary, 1, []);

% Variable tables

V = romsvariablegroups(Opt);

tf = ismember(Opt.vgroups, fieldnames(V));
if ~all(tf)
    str = join(fieldnames(V), ', ');
    error('Unrecognized variable group; options are %s', str);
end

V = cellfun(@(x) V.(x), Opt.vgroups, 'uni', 0);
V = cat(1, V{:});

% Replace all horizontally-dimensioned variables with their boundary
% equivalents

bshort = Opt.boundary;
blong = bshort + "ern boundary condition";

ishvar = ~(endsWith(V.grid, '1') | cellfun(@isempty, V.grid));
Vold = V(ishvar,:);

for ib = length(bshort):-1:1
    Vnew{ib} = Vold;
    Vnew{ib}.short = Vnew{ib}.short + "_" + bshort(ib);
    Vnew{ib}.grid = Vnew{ib}.grid + bshort(ib);
    Vnew{ib}.long = Vnew{ib}.long + ", " + blong(ib);
     
end
V = cat(1, V(~ishvar,:), Vnew{:});

% Initialize schema

Nc = ncschema_init(Opt.format);

% Global attributes

Nc = ncschema_addatts(Nc, 'type', 'BOUNDARY file', ...
                          'history', sprintf('%s: %s', ...
                              datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                              'File schema created via bry_schema.m'));

% Add dimensions

Dims.r1 = {'s_rho'};
Dims.w1 = {'s_w'};

Dims.r2east = {'eta_rho', Opt.tname};
Dims.r3east = {'eta_rho' 's_rho', Opt.tname};
Dims.u2east = {'eta_u', Opt.tname};
Dims.u3east = {'eta_u' 's_rho', Opt.tname};
Dims.v2east = {'eta_v', Opt.tname};
Dims.v3east = {'eta_v' 's_rho', Opt.tname};
Dims.bn3east = {'eta_rho' 'benlayer', Opt.tname};

Dims.r2west = {'eta_rho', Opt.tname};
Dims.r3west = {'eta_rho' 's_rho', Opt.tname};
Dims.u2west = {'eta_u', Opt.tname};
Dims.u3west = {'eta_u' 's_rho', Opt.tname};
Dims.v2west = {'eta_v', Opt.tname};
Dims.v3west = {'eta_v' 's_rho', Opt.tname};
Dims.bn3west = {'eta_rho' 'benlayer', Opt.tname};

Dims.r2south = {'xi_rho', Opt.tname};
Dims.r3south = {'xi_rho' 's_rho', Opt.tname};
Dims.u2south = {'xi_u', Opt.tname};
Dims.u3south = {'xi_u' 's_rho', Opt.tname};
Dims.v2south = {'xi_v', Opt.tname};
Dims.v3south = {'xi_v' 's_rho', Opt.tname};
Dims.bn3south = {'xi_rho' 'benlayer', Opt.tname};

Dims.r2north = {'xi_rho', Opt.tname};
Dims.r3north = {'xi_rho' 's_rho', Opt.tname};
Dims.u2north = {'xi_u', Opt.tname};
Dims.u3north = {'xi_u' 's_rho', Opt.tname};
Dims.v2north = {'xi_v', Opt.tname};
Dims.v3north = {'xi_v' 's_rho', Opt.tname};
Dims.bn3north = {'xi_rho' 'benlayer', Opt.tname};

if Opt.spherical
    coords.r2east = 'lat_rho_east lon_rho_east';
    coords.r3east = 'lat_rho_east lon_rho_east s_rho';
    coords.u2east = 'lat_u_east lon_u_east';
    coords.u3east = 'lat_u_east lon_u_east s_rho';
    coords.v2east = 'lat_v_east lon_v_east';
    coords.v3east = 'lat_v_east lon_v_east s_rho';
    coords.bn3east = 'lat_rho_east lon_rho_east benlayer';

    coords.r2west = 'lat_rho_west lon_rho_west';
    coords.r3west = 'lat_rho_west lon_rho_west s_rho';
    coords.u2west = 'lat_u_west lon_u_west';
    coords.u3west = 'lat_u_west lon_u_west s_rho';
    coords.v2west = 'lat_v_west lon_v_west';
    coords.v3west = 'lat_v_west lon_v_west s_rho';
    coords.bn3west = 'lat_rho_west lon_rho_west benlayer';

    coords.r2south = 'lat_rho_south lon_rho_south';
    coords.r3south = 'lat_rho_south lon_rho_south s_rho';
    coords.u2south = 'lat_u_south lon_u_south';
    coords.u3south = 'lat_u_south lon_u_south s_rho';
    coords.v2south = 'lat_v_south lon_v_south';
    coords.v3south = 'lat_v_south lon_v_south s_rho';
    coords.bn3south = 'lat_rho_south lon_rho_south benlayer';

    coords.r2north = 'lat_rho_north lon_rho_north';
    coords.r3north = 'lat_rho_north lon_rho_north s_rho';
    coords.u2north = 'lat_u_north lon_u_north';
    coords.u3north = 'lat_u_north lon_u_north s_rho';
    coords.v2north = 'lat_v_north lon_v_north';
    coords.v3north = 'lat_v_north lon_v_north s_rho';
    coords.bn3north = 'lat_rho_north lon_rho_north benlayer';
else
    coords.r2east = 'x_rho_east y_rho_east';
    coords.r3east = 'x_rho_east y_rho_east s_rho';
    coords.u2east = 'x_u_east y_u_east';
    coords.u3east = 'x_u_east y_u_east s_rho';
    coords.v2east = 'x_v_east y_v_east';
    coords.v3east = 'x_v_east y_v_east s_rho';
    coords.bn3east = 'x_rho_east y_rho_east benlayer';

    coords.r2west = 'x_rho_west y_rho_west';
    coords.r3west = 'x_rho_west y_rho_west s_rho';
    coords.u2west = 'x_u_west y_u_west';
    coords.u3west = 'x_u_west y_u_west s_rho';
    coords.v2west = 'x_v_west y_v_west';
    coords.v3west = 'x_v_west y_v_west s_rho';
    coords.bn3west = 'x_rho_west y_rho_west benlayer';

    coords.r2south = 'x_rho_south y_rho_south';
    coords.r3south = 'x_rho_south y_rho_south s_rho';
    coords.u2south = 'x_u_south y_u_south';
    coords.u3south = 'x_u_south y_u_south s_rho';
    coords.v2south = 'x_v_south y_v_south';
    coords.v3south = 'x_v_south y_v_south s_rho';
    coords.bn3south = 'x_rho_south y_rho_south benlayer';

    coords.r2north = 'x_rho_north y_rho_north';
    coords.r3north = 'x_rho_north y_rho_north s_rho';
    coords.u2north = 'x_u_north y_u_north';
    coords.u3north = 'x_u_north y_u_north s_rho';
    coords.v2north = 'x_v_north y_v_north';
    coords.v3north = 'x_v_north y_v_north s_rho';
    coords.bn3north = 'x_rho_north y_rho_north benlayer';    
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


