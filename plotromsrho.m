function hp = plotromsrho(Grd, data, axflag)
%PLOTROMSRHO Plots a rho-field variable from ROMS
%
% hp = plotromsrho(Grd, data)
% hp = plotromsrho(Grd, data, flag)
%
% This function is a quick wrapper to create a projected pcolor plot of a
% horizontal slice of ROMS data, using the lat_psi and lon_psi coordinates
% in the grid file structure.   
%
% Input variables:
%
%   Grd:    grid file structure
%
%   data:   nxi x neta array
%
%   flag:   true to issue worldmap command on current axis, false to leave
%           as is (must be map axis)

if nargin < 3
    axflag = true;
end

if axflag
    latlim = minmax(Grd.lat_rho);
    lonlim = minmax(Grd.lon_rho);
    worldmap(latlim, lonlim);
end
    
hp = pcolorm(Grd.lat_psi, Grd.lon_psi, padarray(data(2:end-1,2:end-1), [1 1], NaN, 'post'));
