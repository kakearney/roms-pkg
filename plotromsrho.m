function hp = plotromsrho(Grd, data, axflag)
%PLOTROMSRHO Plots a rho-field variable from ROMS
%
% hp = plotromsrho(Grd, data)
% hp = plotromsrho(Grd, data, flag)
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
