function V = romsvariablegroups(varargin)
%ROMSVARIABLEGROUPS Returns tables of ROMS initialization variables
%
% V = romsvariablegroups;
% V = romsvariablegroups(p1, v1, ...)
%
% This function returns a structure holding tables of ROMS initialization
% variables.  This function is called by various file-schema creation
% functions (ini_schema.m, bry_schema.m) and can also be called directly to
% get info about ROMS I/O variables.  Users can add new categories as
% needed based on ROMS compilation options.
%
% Input variables (optional, passed as parameter/value pairs):
%
%   spherical:  0 or 1, indicating Cartesian or spherical coordinates,
%               respectively [0]
%   
%   Vtransform: 1 or 2, indicating vertical transformation equation (1 =
%               original ROMS transformation, 2 = UCLA flavor
%               transformation) (affects attributes of some reference
%               parameters) [1] 
%
% Output variables:
%
%   V:          structure where each field holds a table with the following
%               columns:
%               short:  character array, variable short name
%               grid:   character array, ROMS grid on which variable is
%                       located:
%                       nul	 non-grided variable
%                       r1   1D rho-variable, s-layer mid-points
%                       w1   1D W-variable, s-layer edges
%                       p2	 2D psi-variable, grid cell corners
%                       r2	 2D rho-variable, grid cell center
%                       u2	 2D U-variable, grid cell x-edges
%                       v2	 2D V-variable, grid cell y-edges
%                       p3	 3D psi-variable, grid cell corners
%                       r3	 3D rho-variable, grid cell center
%                       u3	 3D U-variable, grid cell x-edges
%                       v3	 3D V-variable, grid cell y-edges
%                       w3	 3D W-variable, grid cell center, s-layer edges
%                       b3	 3D BED-sediment variable, grid cell center
%                       l3	 3D spectral light variable, grid cell center
%                       l4	 4D spectral light variable, grid cell center
%                       bn3  3D benthic variable, grid cell center

% Copyright 2022 Kelly Kearney

% Parse inputs

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('spherical', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '<=', 1}));
p.addParameter('Vtransform', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '<=', 2}));
p.addParameter('fillvalue', 1e36, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.parse(varargin{:});
Opt = p.Results;

% Quick helper function to convert cell arrays to tables (more legible to
% write this way than directly constructing the tables)

c2t = @(x) cell2table(x, 'variablenames', {'short', 'grid', 'class', 'long', 'units', 'atts'});

% Variable tables

V.reference = c2t({...
     'spherical'    'nul'  'int32'  'grid type logical switch'                           ''      {'flag_values', [0 1], 'flag_meanings', '0 = Cartesian, 1 = spherical'}
     'Vtransform'   'nul'  'int32'  'vertical terrain-following transformation equation' ''      {}
     'Vstretching'  'nul'  'int32'  'vertical terrain-following stretching function'     ''      {}
     'theta_s'      'nul'  'double' 'S-coordinate surface control parameter'             ''      {}
     'theta_b'      'nul'  'double' 'S-coordinate bottom control parameter'              ''      {}
     'Tcline'       'nul'  'double' 'S-coordinate surface/bottom layer width'            'meter' {}   
     'hc'           'nul'  'double' 'S-coordinate parameter, critical depth'             'meter' {} 
     });


V.grid = c2t({...
     's_rho'        'r1' 'double'  'S-coordinate at RHO-points'                     '' {'valid_min', -1, ...
                                                                                        'valid_max',  0, ...
                                                                                        'positive', 'up', ...
                                                                                        'standard_name', sprintf('ocean_s_coordinate_g%d', Opt.Vtransform), ...
                                                                                        'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'}
     's_w'          'w1' 'double'  'S-coordinate at W-points'                       '' {'valid_min', -1, ...
                                                                                        'valid_max',  0, ...
                                                                                        'positive', 'up', ...
                                                                                        'standard_name', sprintf('ocean_s_coordinate_g%d', Opt.Vtransform), ...
                                                                                        'formula_terms', 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'}
     'Cs_r'         'r1' 'double'  'S-coordinate stretching function at RHO-points' '' {'valid_min', -1, ...
                                                                                        'valid_max',  0}
     'Cs_w'         'w1' 'double'  'S-coordinate stretching function at W-points'   '' {'valid_min', -1, ...
                                                                                        'valid_max',  0}
     });
 
if Opt.spherical
     V.grid = [V.grid; c2t({...
     'h'            'r2' 'double'  'bathymetry at RHO-points'                       'meter' {}
     'lon_rho'      'r2' 'double'  'longitude of RHO-points'                        'degree_east'  {'standard_name', 'longitude'}
     'lat_rho'      'r2' 'double'  'latitude of RHO-points'                         'degree_north' {'standard_name', 'latitude'}
     'lon_u'        'u2' 'double'  'longitude of U-points'                          'degree_east'  {'standard_name', 'longitude'}
     'lat_u'        'u2' 'double'  'latitude of U-points'                           'degree_north' {'standard_name', 'latitude'}
     'lon_v'        'v2' 'double'  'longitude of V-points'                          'degree_east'  {'standard_name', 'longitude'}
     'lat_v'        'v2' 'double'  'latitude of V-points'                           'degree_north' {'standard_name', 'latitude'}
     })];
else
     V.grid = [V.grid; c2t({...
     'h'            'r2' 'double'  'bathymetry at RHO-points'                       'meter' {}
     'x_rho'        'r2' 'double'  'X-location of RHO-points'                       'meter' {}
     'y_rho'        'r2' 'double'  'Y-location of RHO-points'                       'meter' {}
     'x_u'          'u2' 'double'  'X-location of U-points'                         'meter' {}
     'y_u'          'u2' 'double'  'Y-location of U-points'                         'meter' {}
     'x_v'          'v2' 'double'  'X-location of V-points'                         'meter' {}
     'y_v'          'v2' 'double'  'Y-location of V-points'                         'meter' {}
    })];
end
      
V.ocean = c2t({...
     'zeta'         'r2' 'double'  'free-surface'                               'meter'          {'_FillValue', Opt.fillvalue}
     'ubar'         'u2' 'double'  'vertically integrated u-momentum component' 'meter second-1' {'_FillValue', Opt.fillvalue}
     'vbar'         'v2' 'double'  'vertically integrated v-momentum component' 'meter second-1' {'_FillValue', Opt.fillvalue}
     'u'            'u3' 'double'  'u-momentum component'                       'meter second-1' {'_FillValue', Opt.fillvalue}
     'v'            'v3' 'double'  'v-momentum component'                       'meter second-1' {'_FillValue', Opt.fillvalue}
     'temp'         'r3' 'double'  'potential temperature'                      'Celsius'        {'_FillValue', Opt.fillvalue}
     'salt'         'r3' 'double'  'salinity'                                   ''               {'_FillValue', Opt.fillvalue}
     });

V.ice = c2t({...                                                                           
     'uice'         'u2' 'double'  'u-component of ice velocity'                  'meter second-1' {'_FillValue', Opt.fillvalue}
     'vice'         'v2' 'double'  'v-component of ice velocity'                  'meter second-1' {'_FillValue', Opt.fillvalue}
     'aice'         'r2' 'double'  'fraction of cell covered by ice'              ''               {'_FillValue', Opt.fillvalue}
     'hice'         'r2' 'double'  'average ice thickness in cell'                'meter'          {'_FillValue', Opt.fillvalue}
     'snow_thick'   'r2' 'double'  'thickness of snow cover'                      'meter'          {'_FillValue', Opt.fillvalue}
     'tisrf'        'r2' 'double'  'temperature of ice surface'                   'Celsius'        {'_FillValue', Opt.fillvalue}
     'ageice'       'r2' 'double'  'age of the ice'                               'sec'            {'_FillValue', Opt.fillvalue}
     'ti'           'r2' 'double'  'interior ice temperature'                     'Celsius'        {'_FillValue', Opt.fillvalue}
     'sig11'        'r2' 'double'  'internal ice stress 11 component'             'Newton meter-1' {'_FillValue', Opt.fillvalue} 
     'sig22'        'r2' 'double'  'internal ice stress 22 component'             'Newton meter-1' {'_FillValue', Opt.fillvalue} 
     'sig12'        'r2' 'double'  'internal ice stress 12 component'             'Newton meter-1' {'_FillValue', Opt.fillvalue} 
     's0mk'         'r2' 'double'  'salinity of molecular sub-layer under ice'    ''               {'_FillValue', Opt.fillvalue}
     't0mk'         'r2' 'double'  'temperature of molecular sub-layer under ice' 'Celsius'        {'_FillValue', Opt.fillvalue}
     });    
 
 V.BEST_NPZ = c2t({...
     'NO3'          'r3' 'double'  'Nitrate concentration'                        'mmol N m^-3'    {'_FillValue', Opt.fillvalue}
     'NH4'          'r3' 'double'  'Ammonium concentration'                       'mmol N m^-3'    {'_FillValue', Opt.fillvalue}
     'PhS'          'r3' 'double'  'Small phytoplankton concentration'            'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'PhL'          'r3' 'double'  'Large phytoplankton concentration'            'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'MZL'          'r3' 'double'  'Microzooplankton concentration'               'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'Cop'          'r3' 'double'  'Small copepod concentration'                  'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'NCaS'         'r3' 'double'  'On-shelf large copepod concentration'         'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'EupS'         'r3' 'double'  'On-shelf euphausiid concentration'            'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'NCaO'         'r3' 'double'  'Offshore large copepod concentration'         'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'EupO'         'r3' 'double'  'Offshore euphausiid concentration'            'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'Det'          'r3' 'double'  'Slow-sinking detritus concentration'          'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'DetF'         'r3' 'double'  'Fast-sinking detritus concentration'          'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'Jel'          'r3' 'double'  'Jellyfish concentration'                      'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'Ben'          'bn3' 'double' 'Benthic infauna concentration'                'mg C m^-2'      {'_FillValue', Opt.fillvalue}
     'DetBen'       'bn3' 'double' 'Benthic detritus concentration'               'mg C m^-2'      {'_FillValue', Opt.fillvalue}
     'IcePhL'       'r2' 'double'  'Ice algae concentration'                      'mg C m^-3'      {'_FillValue', Opt.fillvalue}
     'IceNO3'       'r2' 'double'  'Ice nitrate concentration'                    'mmol N m^-3'    {'_FillValue', Opt.fillvalue}
     'IceNH4'       'r2' 'double'  'Ice ammonium concentration'                   'mmol N m^-3'    {'_FillValue', Opt.fillvalue}
     'Fe'           'r3' 'double'  'Iron concentration'                           'umol Fe m^-3'   {'_FillValue', Opt.fillvalue}
     'alkalinity'   'r3' 'double'  'total alkalinity'                             'mmol C m^-3'    {'_FillValue', Opt.fillvalue}
     'TIC'          'r3' 'double'  'total inorganic carbon'                       'mmol C m^-3'    {'_FillValue', Opt.fillvalue}
     'oxygen'       'r3' 'double'  'dissolved oxygen concentration'               'mmol O m^-3'    {'_FillValue', Opt.fillvalue}    
    });
 
V.BIO_COBALT = c2t({...
     'nsm'            'r3' 'double' 'Small Phytoplankton Nitrogen'              'mol/kg'           {'_FillValue', Opt.fillvalue}
     'nlg'            'r3' 'double' 'Large Phytoplankton Nitrogen'              'mol/kg'           {'_FillValue', Opt.fillvalue}
     'ndi'            'r3' 'double' 'Diazotroph Nitrogen'                       'mol/kg'           {'_FillValue', Opt.fillvalue}
     'nsmz'           'r3' 'double' 'Small Zooplankton Nitrogen'                'mol/kg'           {'_FillValue', Opt.fillvalue}
     'nmdz'           'r3' 'double' 'Medium-sized zooplankton Nitrogen'         'mol/kg'           {'_FillValue', Opt.fillvalue}
     'nlgz'           'r3' 'double' 'large Zooplankton Nitrogen'                'mol/kg'           {'_FillValue', Opt.fillvalue}
     'ldon'           'r3' 'double' 'labile DON'                                'mol/kg'           {'_FillValue', Opt.fillvalue}
     'sldon'          'r3' 'double' 'Semilabile DON'                            'mol/kg'           {'_FillValue', Opt.fillvalue}
     'srdon'          'r3' 'double' 'Semi-Refractory DON'                       'mol/kg'           {'_FillValue', Opt.fillvalue}
     'nbact'          'r3' 'double' 'bacterial'                                 'mol/kg'           {'_FillValue', Opt.fillvalue}
     'nh4'            'r3' 'double' 'Ammonia'                                   'mol/kg'           {'_FillValue', Opt.fillvalue}
     'no3'            'r3' 'double' 'Nitrate'                                   'mol/kg'           {'_FillValue', Opt.fillvalue}
     'ndet'           'r3' 'double' 'ndet'                                      'mol/kg'           {'_FillValue', Opt.fillvalue}
     'sio4'           'r3' 'double' 'Silicate'                                  'mol/kg'           {'_FillValue', Opt.fillvalue}
     'silg'           'r3' 'double' 'Large Phytoplankton Silicon'               'mol/kg'           {'_FillValue', Opt.fillvalue}
     'sidet'          'r3' 'double' 'Detrital Silicon'                          'mol/kg'           {'_FillValue', Opt.fillvalue}
     'cadet_calc'     'r3' 'double' 'Detrital CaCO3'                            'mol/kg'           {'_FillValue', Opt.fillvalue}
     'cadet_arag'     'r3' 'double' 'Detrital CaCO3'                            'mol/kg'           {'_FillValue', Opt.fillvalue}
     'lith'           'r3' 'double' 'Lithogenic Aluminosilicate'                'g/kg'             {'_FillValue', Opt.fillvalue}
     'lithdet'        'r3' 'double' 'lithdet'                                   'g/kg'             {'_FillValue', Opt.fillvalue}
     'po4'            'r3' 'double' 'Phosphate'                                 'mol/kg'           {'_FillValue', Opt.fillvalue}
     'ldop'           'r3' 'double' 'labile DOP'                                'mol/kg'           {'_FillValue', Opt.fillvalue}
     'sldop'          'r3' 'double' 'Semilabile DOP'                            'mol/kg'           {'_FillValue', Opt.fillvalue}
     'srdop'          'r3' 'double' 'Semi-Refractory DOP'                       'mol/kg'           {'_FillValue', Opt.fillvalue}
     'pdet'           'r3' 'double' 'Detrital Phosphorus'                       'mol/kg'           {'_FillValue', Opt.fillvalue}
     'fesm'           'r3' 'double' 'Small Phytoplankton Iron'                  'mol/kg'           {'_FillValue', Opt.fillvalue}
     'fedi'           'r3' 'double' 'Diazotroph Iron'                           'mol/kg'           {'_FillValue', Opt.fillvalue}
     'felg'           'r3' 'double' 'Large Phytoplankton Iron'                  'mol/kg'           {'_FillValue', Opt.fillvalue}
     'fed'            'r3' 'double' 'Dissolved Iron'                            'mol/kg'           {'_FillValue', Opt.fillvalue}
     'fedet'          'r3' 'double' 'Detrital Iron'                             'mol/kg'           {'_FillValue', Opt.fillvalue}
     'o2'             'r3' 'double' 'Oxygen'                                    'mol/kg'           {'_FillValue', Opt.fillvalue}
     'dic'            'r3' 'double' 'Dissolved Inorganic Carbon'                'mol/kg'           {'_FillValue', Opt.fillvalue}
     'alk'            'r3' 'double' 'Alkalinity'                                'mol/kg'           {'_FillValue', Opt.fillvalue}
     'nmd'            'r3' 'double' 'Medium Phytoplankton Nitrogen'             'mol/kg'           {'_FillValue', Opt.fillvalue}
     'simd'           'r3' 'double' 'Medium Phytoplankton Silicon'              'mol/kg'           {'_FillValue', Opt.fillvalue}
     'femd'           'r3' 'double' 'Medium Phytoplankton Iron'                 'mol/kg'           {'_FillValue', Opt.fillvalue}
     'cased'          'r3' 'double' 'Sediment CaCO3'                            'mol m-3'          {'_FillValue', Opt.fillvalue}
     'cadet_arag_btf' 'r2' 'double' 'aragonite flux to Sediments'               'mol m-2 s-1'      {'_FillValue', Opt.fillvalue}
     'cadet_calc_btf' 'r2' 'double' 'calcite flux to Sediments'                 'mol m-2 s-1'      {'_FillValue', Opt.fillvalue}
     'ndet_btf'       'r2' 'double' 'N flux to Sediments'                       'mol m-2 s-1'      {'_FillValue', Opt.fillvalue}
     'pdet_btf'       'r2' 'double' 'P flux to Sediments'                       'mol m-2 s-1'      {'_FillValue', Opt.fillvalue}
     'sidet_btf'      'r2' 'double' 'SiO2 flux to Sediments'                    'mol m-2 s-1'      {'_FillValue', Opt.fillvalue}
     'chl'            'r3' 'double' 'Chlorophyll'                               'ug/kg'            {'_FillValue', Opt.fillvalue}
     'irr_mem'        'r3' 'double' 'Irradiance memory'                         'Watts/m^2'        {'_FillValue', Opt.fillvalue}
     'htotal'         'r3' 'double' 'H+ ion concentration'                      'mol/kg'           {'_FillValue', Opt.fillvalue}
     'co3_ion'        'r3' 'double' 'Carbonate ion'                             'mol/kg'           {'_FillValue', Opt.fillvalue}
     'mu_mem_sm'      'r3' 'double' 'Aggreg memory small phyto'                 ''                 {'_FillValue', Opt.fillvalue}
     'mu_mem_di'      'r3' 'double' 'Aggreg memory diatoms'                     ''                 {'_FillValue', Opt.fillvalue}
     'mu_mem_lg'      'r3' 'double' 'Aggreg memory large phyto'                 ''                 {'_FillValue', Opt.fillvalue}
     'mu_mem_md'      'r3' 'double' 'Aggreg memory medium phyto'                ''                 {'_FillValue', Opt.fillvalue}
    });

V.BIO_BANAS =  c2t({...
     'phyto'          'r3' 'double'  'phytoplankton'                    'mmol N m^-3' {'_FillValue', Opt.fillvalue}
     'microzoo'       'r3' 'double'  'microzooplankton'                 'mmol N m^-3' {'_FillValue', Opt.fillvalue}
     'det_small'      'r3' 'double'  'small detritus'                   'mmol N m^-3' {'_FillValue', Opt.fillvalue}
     'det_large'      'r3' 'double'  'large detritus (aggregates)'      'mmol N m^-3' {'_FillValue', Opt.fillvalue}
     'NH4'            'r3' 'double'  'ammonium'                         'mmol N m^-3' {'_FillValue', Opt.fillvalue}
     'NO3'            'r3' 'double'  'nitrate'                          'mmol N m^-3' {'_FillValue', Opt.fillvalue}
    });

V.commonbgc = c2t({...
    'ones'            'r3' 'double' 'generic tracer'                   'tracer units' {'_FillValue', Opt.fillvalue}
    'NH4'             'r3' 'double' 'ammonium'                         'umol/kg'      {'_FillValue', Opt.fillvalue}
    'NO3'             'r3' 'double' 'nitrate'                          'umol/kg'      {'_FillValue', Opt.fillvalue}
    'SiO4'            'r3' 'double' 'silicate'                         'umol/kg'      {'_FillValue', Opt.fillvalue}
    'PO4'             'r3' 'double' 'phosphate'                        'umol/kg'      {'_FillValue', Opt.fillvalue}
    'Fe'              'r3' 'double' 'dissolved iron'                   'nmol/L'       {'_FillValue', Opt.fillvalue}
    'TIC'             'r3' 'double' 'total/dissolved inorganic carbon' 'umol/kg'      {'_FillValue', Opt.fillvalue}
    'alkalinity'      'r3' 'double' 'alkalinity'                       'umol/kg'      {'_FillValue', Opt.fillvalue}
    'oxygen'          'r3' 'double' 'oxygen'                           'umol/kg'      {'_FillValue', Opt.fillvalue}
    });


V.nudge = c2t({...
    'M2_NudgeCoef'     'r2' 'double' '2D momentum inverse nudging coefficients'     'day-1' {'_FillValue', Opt.fillvalue}
    'M3_NudgeCoef'     'r3' 'double' '3D momentum inverse nudging coefficients'     'day-1' {'_FillValue', Opt.fillvalue}
    'tracer_NudgeCoef' 'r3' 'double' 'generic tracer inverse nudging coefficients'  'day-1' {'_FillValue', Opt.fillvalue}              
    'temp_NudgeCoef'   'r3' 'double' 'temperature inverse nudging coefficients'     'day-1' {'_FillValue', Opt.fillvalue}
    'salt_NudgeCoef'   'r3' 'double' 'salinity inverse nudging coefficients'        'day-1' {'_FillValue', Opt.fillvalue}
    });

% Fix the attribute formatting, where necessary

V = structfun(@fixatts, V, 'uni', 0);


function tbl = fixatts(tbl)
%FIXATTS Fix attribute formatting
%
% Enforce atts column to be an nx1 cell array, even if all rows include
% same-sized cell arrays.

if size(tbl.atts,2)>1
    for ii = 1:height(tbl)
        tbl.atts{ii,1} =  tbl.atts(ii,:);
    end
    tbl.atts = tbl.atts(:,1);
end