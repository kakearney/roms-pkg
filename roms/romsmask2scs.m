function Scs = romsmask2scs(rhomask)
%ROMSMASK2SCS Returns start/count/stride structure for a mask
%
% This function returns the start/count/stride data defining the horizontal
% hyperslab emcompassing a rho-grid mask for a ROMS grid.  This includes
% the extent of true-masked rho-points plus a one-grid-cell buffer, along
% with the psi-, u-, and v-points encompassed by that mask.  The intent is
% to allow easy reading of a sub-grid of cells with the same
% characteristics as a full roms grid, to allow for easy analysis and
% plotting by functions expecting input to conform to those size
% conventions.
%
% Input variables:
%
%   rhomask: nxi x neta logical array, where nxi and neta correspond to the
%            size of the ROMS rho grid
%
% Output variables:
%
%   Scs:    1 x 1 structure with field names corresponding to ROMS grid
%           horizontal grid dimensions.  Each hold [start count stride]
%           information that can be passed to netCDF-reading functions
%           (see ncstruct, ncread, etc.) to read a specific hyperslab.

% Copyright 2024 Kelly Kearney

[xx,ee] = ind2sub(size(rhomask), find(rhomask));
xlim = minmax(xx);
elim = minmax(ee);

nx = diff(xlim)+1;
ne = diff(elim)+1;

Scs = struct('xi_rho',  [xlim(1)-1 nx+2 1], ...
             'eta_rho', [elim(1)-1 ne+2 1], ...
             'xi_psi',  [xlim(1)-1 nx+1 1], ...
             'eta_psi', [elim(1)-1 ne+1 1], ...
             'xi_u',    [xlim(1)-1 nx+1 1], ...
             'eta_u',   [elim(1)-1 nx+2 1], ...
             'xi_v',    [xlim(1)-1 nx+2 1], ...
             'eta_v',   [elim(1)-1 ne+1 1]);
