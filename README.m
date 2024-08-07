%% Some ROMS utilities
% Author: Kelly Kearney
%
% This repository holds a loose collections of utilities I've written to
% interact with output from the <https://www.myroms.org Regional Ocean
% Modeling System (ROMS)>, a widely-used primitive-equations ocean model.
%
%% Getting started
%
% *Prerequisites*
%
% For full functionality, these tools require Matlab R2012a or later with
% the Mapping, Statistics, and Image Processing Toolboxes.
%
% *Downloading and installation*
%
% This code can be downloaded from <https://github.com/kakearney/roms-pkg/ Github>. 
%
% *Matlab Search Path*
%
% The following folders need to be added to your Matlab Search path (via
% |addpath|, |pathtool|, etc.):
%
%  roms-pkg/roms
%  roms-pkg/cellstr2
%  roms-pkg/invdisweight
%  roms-pkg/mask2poly
%  roms-pkg/minmax
%  roms-pkg/regexpfound
%  roms-pkg/setgetpos_V1.2
%  roms-pkg/ConsoleProgressBar
%
% The primary functions for this toolbox can be found in the roms folder; 
% the remaining folders hold dependent functions used by thse primary ones.
%
% In addition to the dependencies listed above, some roms-pkg functions rely
% on the netCDF-reading utilities from the
% <https://github.com/chadagreene/CDT Climate Data Toolbox>; I 
% highly recommend users download this entire toolbox for its extra
% netCDF utilities and its variety of plotting and analysis functions that
% can aid in working with ROMS output datasets.  Finally, the ROMS folder
% includes two utility functions, |stretching| and |set_depth|, that are
% copied over from the <https://www.myroms.org/wiki/Matlab_Scripts
% myroms.org Matlab Scripts>; personally, I'm not a fan of this toolbox due
% to its clunky syntax, but it's useful to have handy occasionally if you
% work with ROMS output.

%% Functions
%
% * |plotromsrho|: Plot a 2D slice of a rho variable to projected pcolor
%   plot
% * |animateromsrho|: Animate a 2D slice of a rho variable over time
% * |calcromsz|: Calculate rho- and w-depth values based on bottom depth
%   and surface height
% * |romsboundarycoords|: Extract boundary slices from 4D grid
% * |romsgeometryparams|: Calculate various grid-related geometry
%   parameters
% * |parsevarinfo|: Read data from a varinfo.dat file into a table
% * |parseromslog|: Read archived standard output text from from a ROMS
%   simulation
% * |oceandata2romsgrd|: Regrid global gridded ocean data to a ROMS grid
% * |bry_schema|: Build a netCDF file schema for a ROMS boundary file
% * |ini_schema|: Build a netCDF file schema for a ROMS initialization file
% * |nud_schema|: Build a netCDF file schema for a ROMS nudging file
% * |romsavgclimatology|: Calculate climatology from ROMS averages output
% * |romsmask2scs|: Returns netCDF start/count/stride hyperslab info
%   (across all ROMS dimensions) encompassing all true values in a rho-grid
%   mask
% * |romstransect|: extracts xi- and eta-vs-depth transects from 3D ROMS
%   output in distance-along-transect (km) and depth (m) coordinates
% * |romsvariablegroups|: Returns tables of ROMS variables associated with
%   common modules

%% Contributions
%
% Community contributions to this package are welcome!
% 
% To report bugs, please submit
% <https://github.com/kakearney/example-pkg/issues an issue> on GitHub and
% include:  
% 
% * your operating system
% * your version of Matlab and all relevant toolboxes (type |ver| at the Matlab command line to get this info)  
% * code/data to reproduce the error or buggy behavior, and the full text of any error messages received 
% 
% Please also feel free to submit enhancement requests, or to send pull
% requests (via GitHub) for bug fixes or new features. 
% 
% I do monitor the MatlabCentral FileExchange entry for any issues raised
% in the comments, but would prefer to track issues on GitHub. 
% 

