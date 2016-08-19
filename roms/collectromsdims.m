function Dim = collectromsdims(files, type, nz, h, flag)
%COLLECTROMSDIMS Read time and depth data from ROMS his, avg, or sta files
%
% Dim = collectromsdims(files, type, nz, h)
% Dim = collectromsdims(files, type, nz, h, flag)
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
%   h:      nxi x neta array of bottom depth values (his/avg only)
%
%   flag:   logical scalar, true to calculate zr and zw (can be
%           time-consuming for large runs).  Default is true.
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

if nargin < 5
    flag = true;
end

switch type
    
    case {'his', 'avg'}
        if ischar(files)
            files = {files};
        end
        
        S = ncreads(files{1}, 'Vstretching', 'Vtransform', 'theta_s', 'theta_b', 'hc');
        
        zeta = cellfun(@(x) ncread(x, 'zeta'), files, 'uni', 0);
        Dim.zeta = cat(3, zeta{:});
        
        t = cellfun(@(x) ncread(x, 'ocean_time'), files, 'uni', 0);
        Dim.ocean_time = cat(1, t{:});
        
        tunit = ncreadatt(files{1}, 'ocean_time', 'units');
        tparts = textscan(tunit, '%s since %D', 1);
        switch lower(tparts{1}{1})
            case 'seconds'
                Dim.date = tparts{2} + seconds(Dim.ocean_time);
            case 'hours'
                Dim.date = tparts{2} + hours(Dim.ocean_time);
            case 'days'
                Dim.date = tparts{2} + days(Dim.ocean_time);
            otherwise
                warning('Could not parse reference time');
                Dim.date = [];
        end
        
        if flag
            [zr, zw] = cellfun(@(x) calcromsz(h, x, nz, ...
                'Vstretching', S.Vstretching, ...
                'Vtransform', S.Vtransform, ...
                'theta_s', S.theta_s, ...
                'theta_b', S.theta_b, ...
                'hc', S.hc), zeta, 'uni', 0);
            Dim.zr = cat(4, zr{:});
            Dim.zw = cat(4, zw{:});
            
%             tic
%             [zr, zw] = calcromsz(h, Dim.zeta, nz, ...
%                 'Vstretching', S.Vstretching, ...
%                 'Vtransform', S.Vtransform, ...
%                 'theta_s', S.theta_s, ...
%                 'theta_b', S.theta_b, ...
%                 'hc', S.hc);
%             toc

        else
            Dim.zr = [];
            Dim.zw = [];
        end
        
%         Dim.zeta = zeta;
%         Dim.zr = zr;
%         Dim.zw = zw;
        
    case 'sta'
        
        if ischar(files)
           files = {files}; 
        end
        
        S = ncreads(files{1}, 'Vstretching', 'Vtransform', 'theta_s', 'theta_b', 'hc');
        
        Dim = ncreads(files{1}, 'lat_rho', 'lon_rho', 'h', 'ocean_time', 'zeta');
        
        Tmp = cellfun(@(x) ncreads(x, 'ocean_time', 'zeta'), files);
        
        Dim.ocean_time = cat(1, Tmp.ocean_time);
        Dim.zeta = cat(2, Dim.zeta);
        
        
        tunit = ncreadatt(files{1}, 'ocean_time', 'units');
        tparts = textscan(tunit, '%s since %D', 1);
        switch lower(tparts{1}{1})
            case 'seconds'
                Dim.date = tparts{2} + seconds(Dim.ocean_time);
            case 'hours'
                Dim.date = tparts{2} + hours(Dim.ocean_time);
            case 'days'
                Dim.date = tparts{2} + days(Dim.ocean_time);
            otherwise
                warning('Could not parse reference time');
                Dim.date = [];
        end
        
        
        if flag
            [zr, zw] = cellfun(@(x) calcromsz(Dim.h, x, nz, ...
                'Vstretching', S.Vstretching, ...
                'Vtransform', S.Vtransform, ...
                'theta_s', S.theta_s, ...
                'theta_b', S.theta_b, ...
                'hc', S.hc), {Tmp.zeta}, 'uni', 0);
            
            Dim.zr = permute(cat(4, zr{:}), [3 1 4 2]);
            Dim.zw = permute(cat(4, zw{:}), [3 1 4 2]);
        else
            Dim.zr = [];
            Dim.zw = [];
        end
        
end