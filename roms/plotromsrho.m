function hp = plotromsrho(Grd, data)

% h = plotgrid('setup', cell(1), [],[], 'mar', 0.02);
% h.fig.Color = 'w';

latlim = minmax(Grd.lat_rho);
lonlim = minmax(Grd.lon_rho);

worldmap(latlim, lonlim);
hp = pcolorm(Grd.lat_psi, Grd.lon_psi, padarray(data(2:end-1,2:end-1), [1 1], NaN, 'post'));
% setm(ax, 'flinewidth', 1);