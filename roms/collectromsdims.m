function Dim = collectromsdims(files, type, nz, h)
%COLLECTROMSDIMS Read time and depth data from ROMS his, avg, or sta files
%
% This function read the time and depth data from ROMS history, average, or
% stations output files, and converts to depth space.
%
% Input variables
%
%   files:  file name(s) as string or cell array of strings
%
%   type:   'his', 'avg', or 'sta'
%
%   nz:     number of depth levels in model domain
%
%   h:      nxi x neta array of bottom depth values
%
% Output variables:
%
%   Dim:    1 x 1 structure with the following fields:
%
%           ocean_time: nt x 1 array, model time
%
%           zeta:       nxi x neta x nt array, free surface height
%
%           zr:         nxi x neta x nz x nt array, depth at rho points
%
%           zw:         nxi x neta x nz+1 x nt array, depth at psi points

% Copyright 2016 Kelly Kearney

switch type
    
    case {'his', 'avg'}
        if ischar(files)
            files = {files};
        end
        
        S = ncreads(files{1}, 'Vstretching', 'Vtransform', 'theta_s', 'theta_b', 'hc');
        
        zeta = cellfun(@(x) ncread(x, 'zeta'), files, 'uni', 0);
        zeta = cat(3, zeta{:});
        
        t = cellfun(@(x) ncread(x, 'ocean_time'), files, 'uni', 0);
        Dim.ocean_time = cat(1, t{:});
        
        [zr, zw] = calcromsz(h, zeta, nz, ...
            'Vstretching', S.Vstretching, ...
            'Vtransform', S.Vtransform, ...
            'theta_s', S.theta_s, ...
            'theta_b', S.theta_b, ...
            'hc', S.hc);
        
        Dim.zeta = zeta;
        Dim.zr = zr;
        Dim.zw = zw;
        
    case 'sta'
        
        if ~ischar(files)
            error('Only a single station file can be read at a time');
        end
        
        S = ncreads(files, 'Vstretching', 'Vtransform', 'theta_s', 'theta_b', 'hc');
        Dim = ncreads(files, 'lat_rho', 'lon_rho', 'h', 'ocean_time', 'zeta');
        
        [zr, zw] = calcromsz(Dim.h, Dim.zeta, nz, ...
            'Vstretching', S.Vstretching, ...
            'Vtransform', S.Vtransform, ...
            'theta_s', S.theta_s, ...
            'theta_b', S.theta_b, ...
            'hc', S.hc);
        
        Dim.zr = permute(zr, [3 1 4 2]);
        Dim.zw = permute(zw, [3 1 4 2]);
        
end