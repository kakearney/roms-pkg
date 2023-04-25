function Nc = nud_schema(varargin)
%NUD_SCHEMA Netcdf file schema for a ROMS nudging file
%
%

p = inputParser;
p.addParameter('spherical', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '<=', 1}));
p.addParameter('Vtransform', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '<=', 2}));
p.addParameter('Lp', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('Mp', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('N',  1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('NT', 2, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('NBL', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
p.addParameter('format', 'classic', @(x) validateattributes(x, {'string', 'char'}, {'scalartext'}));
p.addParameter('vgroups', {'reference', 'grid', 'nudge'}, @(x) iscellstr(x) || isstring(x));

p.parse(varargin{:});
Opt = p.Results;

% Variable tables

if ismember('ocean', Opt.vgroups)
    warning('Replacing "ocean" group with "nudge" group');
    Opt.vgroups = union(setdiff(Opt.vgroups, 'ocean'), 'nudge');
end

V = romsvariablegroups(Opt);
 
% Parse variable groups

tf = ismember(Opt.vgroups, fieldnames(V));
if ~all(tf)
    str = join(fieldnames(V), ', ');
    error('Unrecognized variable group; options are %s', str);
end

if ismember('reference', Opt.vgroups)
    V.reference = V.reference(strcmp(V.reference.short, 'spherical'),:);
end
if ismember('grid', Opt.vgroups)
    V.grid = V.grid(contains(V.grid.short, 'rho'),:);
end

extragrp = setdiff(Opt.vgroups, {'reference', 'grid', 'nudge'});
if ~isempty(extragrp)
    Vex = cellfun(@(x) V.(x), extragrp, 'uni', 0);
    Vex = cat(1, Vex{:});

    Vex.short = Vex.short + "_NudgCoef";
    Vex.long = Vex.long + " inverse nudging coefficients";
    Vex.units = repmat({'day-1'}, height(Vex), 1);
    
    V = cellfun(@(x) V.(x), setdiff(Opt.vgroups, extragrp), 'uni', 0);
    V = [cat(1, V{:}); Vex];
else
    V = cellfun(@(x) V.(x), setdiff(Opt.vgroups, extragrp), 'uni', 0);
    V = cat(1, V{:});
end
    
% Initialize schema

Nc = ncschema_init(Opt.format);

% Global attributes

Nc = ncschema_addatts(Nc, 'type', 'Nudging Coefficients file', ...
                          'history', sprintf('%s: %s', ...
                              datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                              'File schema created via nud_schema.m'));

% Add dimensions

Dims.r1 = {'s_rho'};
Dims.w1 = {'s_w'};
Dims.r2 = {'xi_rho' 'eta_rho'};
Dims.r3 = {'xi_rho' 'eta_rho' 's_rho'};
Dims.u2 = {'xi_u' 'eta_u'};
Dims.u3 = {'xi_u' 'eta_u' 's_rho'};
Dims.v2 = {'xi_v' 'eta_v'};
Dims.v3 = {'xi_v' 'eta_v' 's_rho'};
Dims.bn3 = {'xi_rho' 'eta_rho' 'benlayer'};
    
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
    'benlayer',     Opt.NBL,  false);

dused = cellfun(@(x) Dims.(x), setdiff(unique(V.grid), 'nul'), 'uni', 0);
dused = unique(cat(2, dused{:}));

tf = ismember({Nc.Dimensions.Name}, dused);
Nc.Dimensions = Nc.Dimensions(tf);

% Add variables

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
    
    atts = [atts V.atts{iv}];
    
    % Add to schema
    
    if strcmp(V.grid{iv}, 'nul')
        dims = {};
    else
        dims = Dims.(V.grid{iv});
    end
    
    Nc = ncschema_addvars(Nc, V.short{iv}, dims, atts, V.class{iv});
    
end

