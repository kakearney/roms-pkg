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
     'zeta'         'r2' 'double'  'free-surface'                               'meter'          {}
     'ubar'         'u2' 'double'  'vertically integrated u-momentum component' 'meter second-1' {}
     'vbar'         'v2' 'double'  'vertically integrated v-momentum component' 'meter second-1' {}
     'u'            'u3' 'double'  'u-momentum component'                       'meter second-1' {}
     'v'            'v3' 'double'  'v-momentum component'                       'meter second-1' {}
     'temp'         'r3' 'double'  'potential temperature'                      'Celsius'        {}
     'salt'         'r3' 'double'  'salinity'                                   ''               {}
     });

V.ice = c2t({...                                                                           
     'uice'         'u2' 'double'  'u-component of ice velocity'                  'meter second-1' {}
     'vice'         'v2' 'double'  'v-component of ice velocity'                  'meter second-1' {}
     'aice'         'r2' 'double'  'fraction of cell covered by ice'              ''               {}
     'hice'         'r2' 'double'  'average ice thickness in cell'                'meter'          {}
     'snow_thick'   'r2' 'double'  'thickness of snow cover'                      'meter'          {}
     'tisrf'        'r2' 'double'  'temperature of ice surface'                   'Celsius'        {}
     'ageice'       'r2' 'double'  'age of the ice'                               'sec'            {}
     'ti'           'r2' 'double'  'interior ice temperature'                     'Celsius'        {}
     'sig11'        'r2' 'double'  'internal ice stress 11 component'             'Newton meter-1' {} 
     'sig22'        'r2' 'double'  'internal ice stress 22 component'             'Newton meter-1' {} 
     'sig12'        'r2' 'double'  'internal ice stress 12 component'             'Newton meter-1' {} 
     's0mk'         'r2' 'double'  'salinity of molecular sub-layer under ice'    ''               {}
     't0mk'         'r2' 'double'  'temperature of molecular sub-layer under ice' 'Celsius'        {}
     });
    
 
 V.BEST_NPZ = c2t({...
     'NO3'          'r3' 'double'  'Nitrate concentration'                        'mmol N m^-3'    {}
     'NH4'          'r3' 'double'  'Ammonium concentration'                       'mmol N m^-3'    {}
     'PhS'          'r3' 'double'  'Small phytoplankton concentration'            'mg C m^-3'      {}
     'PhL'          'r3' 'double'  'Large phytoplankton concentration'            'mg C m^-3'      {}
     'MZL'          'r3' 'double'  'Microzooplankton concentration'               'mg C m^-3'      {}
     'Cop'          'r3' 'double'  'Small copepod concentration'                  'mg C m^-3'      {}
     'NCaS'         'r3' 'double'  'On-shelf large copepod concentration'         'mg C m^-3'      {}
     'EupS'         'r3' 'double'  'On-shelf euphausiid concentration'            'mg C m^-3'      {}
     'NCaO'         'r3' 'double'  'Offshore large copepod concentration'         'mg C m^-3'      {}
     'EupO'         'r3' 'double'  'Offshore euphausiid concentration'            'mg C m^-3'      {}
     'Det'          'r3' 'double'  'Slow-sinking detritus concentration'          'mg C m^-3'      {}
     'DetF'         'r3' 'double'  'Fast-sinking detritus concentration'          'mg C m^-3'      {}
     'Jel'          'r3' 'double'  'Jellyfish concentration'                      'mg C m^-3'      {}
     'Ben'          'bn3' 'double'  'Benthic infauna concentration'                'mg C m^-2'      {}
     'DetBen'       'bn3' 'double'  'Benthic detritus concentration'               'mg C m^-2'      {}
     'IcePhL'       'r2' 'double'  'Ice algae concentration'                      'mg C m^-3'      {}
     'IceNO3'       'r2' 'double'  'Ice nitrate concentration'                    'mmol N m^-3'    {}
     'IceNH4'       'r2' 'double'  'Ice ammonium concentration'                   'mmol N m^-3'    {}
     'Fe'           'r3' 'double'  'Iron concentration'                           'umol Fe m^-3'   {}
     'alkalinity'   'r3' 'double'  'total alkalinity'                             'mmol C m^-3'    {}
     'TIC'          'r3' 'double'  'total inorganic carbon'                       'mmol C m^-3'    {}
     'oxygen'       'r3' 'double'  'dissolved oxygen concentration'               'mmol O m^-3'    {}    
    });
 
V.BIO_COBALT = c2t({...
     'nsm'            'r3' 'double' 'Small Phytoplankton Nitrogen'              'mol/kg'           {}        
     'nlg'            'r3' 'double' 'Large Phytoplankton Nitrogen'              'mol/kg'           {}         
     'ndi'            'r3' 'double' 'Diazotroph Nitrogen'                       'mol/kg'           {}         
     'nsmz'           'r3' 'double' 'Small Zooplankton Nitrogen'                'mol/kg'           {}         
     'nmdz'           'r3' 'double' 'Medium-sized zooplankton Nitrogen'         'mol/kg'           {}         
     'nlgz'           'r3' 'double' 'large Zooplankton Nitrogen'                'mol/kg'           {}         
     'ldon'           'r3' 'double' 'labile DON'                                'mol/kg'           {}         
     'sldon'          'r3' 'double' 'Semilabile DON'                            'mol/kg'           {}         
     'srdon'          'r3' 'double' 'Semi-Refractory DON'                       'mol/kg'           {}         
     'nbact'          'r3' 'double' 'bacterial'                                 'mol/kg'           {}         
     'nh4'            'r3' 'double' 'Ammonia'                                   'mol/kg'           {}         
     'no3'            'r3' 'double' 'Nitrate'                                   'mol/kg'           {}         
     'ndet'           'r3' 'double' 'ndet'                                      'mol/kg'           {}         
     'sio4'           'r3' 'double' 'Silicate'                                  'mol/kg'           {}         
     'silg'           'r3' 'double' 'Large Phytoplankton Silicon'               'mol/kg'           {}         
     'sidet'          'r3' 'double' 'Detrital Silicon'                          'mol/kg'           {}         
     'cadet_calc'     'r3' 'double' 'Detrital CaCO3'                            'mol/kg'           {}         
     'cadet_arag'     'r3' 'double' 'Detrital CaCO3'                            'mol/kg'           {}         
     'lith'           'r3' 'double' 'Lithogenic Aluminosilicate'                'g/kg'             {}       
     'lithdet'        'r3' 'double' 'lithdet'                                   'g/kg'             {}       
     'po4'            'r3' 'double' 'Phosphate'                                 'mol/kg'           {}         
     'ldop'           'r3' 'double' 'labile DOP'                                'mol/kg'           {}         
     'sldop'          'r3' 'double' 'Semilabile DOP'                            'mol/kg'           {}         
     'srdop'          'r3' 'double' 'Semi-Refractory DOP'                       'mol/kg'           {}         
     'pdet'           'r3' 'double' 'Detrital Phosphorus'                       'mol/kg'           {}         
     'fesm'           'r3' 'double' 'Small Phytoplankton Iron'                  'mol/kg'           {}         
     'fedi'           'r3' 'double' 'Diazotroph Iron'                           'mol/kg'           {}         
     'felg'           'r3' 'double' 'Large Phytoplankton Iron'                  'mol/kg'           {}         
     'fed'            'r3' 'double' 'Dissolved Iron'                            'mol/kg'           {}         
     'fedet'          'r3' 'double' 'Detrital Iron'                             'mol/kg'           {}         
     'o2'             'r3' 'double' 'Oxygen'                                    'mol/kg'           {}         
     'dic'            'r3' 'double' 'Dissolved Inorganic Carbon'                'mol/kg'           {}         
     'alk'            'r3' 'double' 'Alkalinity'                                'mol/kg'           {}         
     'nmd'            'r3' 'double' 'Medium Phytoplankton Nitrogen'             'mol/kg'           {}         
     'simd'           'r3' 'double' 'Medium Phytoplankton Silicon'              'mol/kg'           {}         
     'femd'           'r3' 'double' 'Medium Phytoplankton Iron'                 'mol/kg'           {}         
     'cased'          'r3' 'double' 'Sediment CaCO3'                            'mol m-3'          {}          
     'cadet_arag_btf' 'r2' 'double' 'aragonite flux to Sediments'               'mol m-2 s-1'      {}              
     'cadet_calc_btf' 'r2' 'double' 'calcite flux to Sediments'                 'mol m-2 s-1'      {}              
     'ndet_btf'       'r2' 'double' 'N flux to Sediments'                       'mol m-2 s-1'      {}              
     'pdet_btf'       'r2' 'double' 'P flux to Sediments'                       'mol m-2 s-1'      {}              
     'sidet_btf'      'r2' 'double' 'SiO2 flux to Sediments'                    'mol m-2 s-1'      {}              
     'chl'            'r3' 'double' 'Chlorophyll'                               'ug/kg'            {}        
     'irr_mem'        'r3' 'double' 'Irradiance memory'                         'Watts/m^2'        {}            
     'htotal'         'r3' 'double' 'H+ ion concentration'                      'mol/kg'           {}         
     'co3_ion'        'r3' 'double' 'Carbonate ion'                             'mol/kg'           {}         
     'mu_mem_sm'      'r3' 'double' 'Aggreg memory small phyto'                 ''                 {}   
     'mu_mem_di'      'r3' 'double' 'Aggreg memory diatoms'                     ''                 {}   
     'mu_mem_lg'      'r3' 'double' 'Aggreg memory large phyto'                 ''                 {}   
     'mu_mem_md'      'r3' 'double' 'Aggreg memory medium phyto'                ''                 {}  
    });
 

V.BIO_BANAS =  c2t({...
     'phyto'          'r3' 'double'  'phytoplankton'                    'mmol N m^-3' {}  
     'microzoo'       'r3' 'double'  'microzooplankton'                 'mmol N m^-3' {}
     'det_small'      'r3' 'double'  'small detritus'                   'mmol N m^-3' {}
     'det_large'      'r3' 'double'  'large detritus (aggregates)'      'mmol N m^-3' {}
     'NH4'            'r3' 'double'  'ammonium'                         'mmol N m^-3' {}
     'NO3'            'r3' 'double'  'nitrate'                          'mmol N m^-3' {}
    });

V.commonbgc = c2t({...
    'ones'            'r3' 'double' 'generic tracer'                   'tracer units' {}
    'NH4'             'r3' 'double' 'ammonium'                         'umol/kg' {}
    'NO3'             'r3' 'double' 'nitrate'                          'umol/kg' {}
    'SiO4'            'r3' 'double' 'silicate'                         'umol/kg' {}
    'PO4'             'r3' 'double' 'phosphate'                        'umol/kg' {}
    'Fe'              'r3' 'double' 'dissolved iron'                   'nmol/L'  {}
    'TIC'             'r3' 'double' 'total/dissolved inorganic carbon' 'umol/kg' {}
    'alkalinity'      'r3' 'double' 'alkalinity'                       'umol/kg' {}
    'oxygen'          'r3' 'double' 'oxygen'                           'umol/kg' {}
    });