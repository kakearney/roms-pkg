function B = romsboundarycoords(Grd, S)
%ROMSBOUNDARYCOORDS Calculates coordinates of ROMS boundary conditions
%
% B = romsboundarycoords(Grd, S)
%
% Quick helper function to extract the grid slices associated with each
% boundary condition
%
% Input variables:
%
%   Grd:    structure array corresponding to ROMS grid file (i.e.
%           ncreads(grdfile)
%
%   S:      structure array holding geometry variables (see
%           romsgeometryparams.m)
%
% Output variables:
%
%   B:      structure array whose fields correspond to the latitude,
%           longitude, and depth (m) coordinates of each boundary condition
%           slice.  Fields are names as [var]_[grid]_[bound], where [var]
%           is 'lat', 'lon', or 'z'; [grid] is 'rho', 'u', or 'v'; and
%           [bound] is 'west', 'east', 'south', or 'north'.

% Copyright 2017 Kelly Kearney


% rho points

B.lat_rho_west  = Grd.lat_rho(1  , :  );
B.lat_rho_east  = Grd.lat_rho(end, :  );
B.lat_rho_south = Grd.lat_rho(:  , 1  );
B.lat_rho_north = Grd.lat_rho(:  , end);

B.lon_rho_west  = Grd.lon_rho(1  , :  );
B.lon_rho_east  = Grd.lon_rho(end, :  );
B.lon_rho_south = Grd.lon_rho(:  , 1  );
B.lon_rho_north = Grd.lon_rho(:  , end);

B.z_rho_west  = permute(S.z_r(1,  :,  :), [2 3 1]);
B.z_rho_east  = permute(S.z_r(end,:,  :), [2 3 1]);
B.z_rho_south = permute(S.z_r(:,  1,  :), [1 3 2]);
B.z_rho_north = permute(S.z_r(:,  end,:), [1 3 2]);

% u points

B.lat_u_west  = Grd.lat_u(1  , :  );
B.lat_u_east  = Grd.lat_u(end, :  );
B.lat_u_south = Grd.lat_u(:  , 1  );
B.lat_u_north = Grd.lat_u(:  , end);

B.lon_u_west  = Grd.lon_u(1  , :  );
B.lon_u_east  = Grd.lon_u(end, :  );
B.lon_u_south = Grd.lon_u(:  , 1  );
B.lon_u_north = Grd.lon_u(:  , end);

B.z_u_west  = permute(S.z_u(1,  :,  :), [2 3 1]);
B.z_u_east  = permute(S.z_u(end,:,  :), [2 3 1]);
B.z_u_south = permute(S.z_u(:,  1,  :), [1 3 2]);
B.z_u_north = permute(S.z_u(:,  end,:), [1 3 2]);

% v points

B.lat_v_west  = Grd.lat_v(1  , :  );
B.lat_v_east  = Grd.lat_v(end, :  );
B.lat_v_south = Grd.lat_v(:  , 1  );
B.lat_v_north = Grd.lat_v(:  , end);

B.lon_v_west  = Grd.lon_v(1  , :  );
B.lon_v_east  = Grd.lon_v(end, :  );
B.lon_v_south = Grd.lon_v(:  , 1  );
B.lon_v_north = Grd.lon_v(:  , end);

B.z_v_west  = permute(S.z_v(1,  :,  :), [2 3 1]);
B.z_v_east  = permute(S.z_v(end,:,  :), [2 3 1]);
B.z_v_south = permute(S.z_v(:,  1,  :), [1 3 2]);
B.z_v_north = permute(S.z_v(:,  end,:), [1 3 2]);







