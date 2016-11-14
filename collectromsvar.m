function [bhis, bsta] = collectromsvar(hisfiles, stafiles, var, iz)
%COLLECTROMSVAR Read data from history and station files
%
% [bhis, bsta] = collectromsvar(hisfiles, stafiles, var, iz)
%
% This function provides a way to read a variable from the potentially
% numerous history and station files that make up a ROMS simulation,
% without needing to keep track of dimensions, 3D vs 4D variables, etc.
%
% Input variables:
%
%   hisfiles:   string or cell array of strings, names of history files
%
%   stafiles:   string or cell array of strings, names of station files
%
%   var:        name of variable to extract
%
%   iz:         depth layers to extract.  Should be sorted.
%
% Output variables:
%
%   bhis:       nxi x neta x nt (x nz) array of data from history files
%
%   bsta:       nstat x nt (x nz) array of data from station files

% Copyright 2016 Kelly Kearney

% Check that data files exist

if ischar(hisfiles)
    hisfiles = {hisfiles};
end
if ischar(stafiles)
    stafiles = {stafiles};
end
hisfiles = hisfiles(:);
stafiles = stafiles(:);

fexists = cellfun(@(x) exist(x,'file'), [hisfiles; stafiles]);
if ~all(fexists)
    tmp = [hisfiles; stafiles];
    tmp = sprintf('  %s\n', tmp{~fexists});
    error('File(s) not found:\n%s', tmp);
end

% Read some file info

if isempty(hisfiles)
    Tmp = ncinfo(stafiles{1});
    isz = strcmp({Tmp.Dimensions.Name}, 's_rho');
    nz = Tmp.Dimensions(isz).Length;
     
    isvar = strcmp({Tmp.Variables.Name}, var);
    ndim = length(Tmp.Variables(isvar).Size)+1;
    
else
    Tmp = ncinfo(hisfiles{1});
    isz = strcmp({Tmp.Dimensions.Name}, 'N');
    if ~any(isz)
        isz = strcmp({Tmp.Dimensions.Name}, 's_rho');
    end
    nz = Tmp.Dimensions(isz).Length;

    isvar = strcmp({Tmp.Variables.Name}, var);
    if ~any(isvar)
        error('Variable (%s) not found in file', var);
    end
    ndim = length(Tmp.Variables(isvar).Size);
end

% Parse start/count/step from iz, to read minimum possible

if ndim == 3
    htmp = cellfun(@(x) ncread(x, var), hisfiles, 'uni', 0);
    stmp = cellfun(@(x) ncread(x, var), stafiles, 'uni', 0);
    
    bhis = cat(3, htmp{:});
    bsta = cat(2, stmp{:});
    
    bsta(bsta > 1e35) = NaN;
    
elseif ndim == 4
    
    flag = false;
    if isscalar(iz)
        stt = iz;
        cnt = 1;
        stp = 1;
    else
        iz = sort(iz);
        dlayer = diff(iz);
        if length(unique(dlayer)) == 1
            stt = min(iz);
            cnt = length(iz);
            stp = dlayer(1);
        else % Some assortment, read all and then extract
            stt = 1;
            cnt = Inf;
            stp = 1;
            flag = true;
        end
    end
    
    htmp = cellfun(@(x) ncread(x, var, [1 1 stt 1], [Inf Inf cnt Inf], [1 1 stp 1]), hisfiles, 'uni', 0);
    stmp = cellfun(@(x) ncread(x, var, [stt 1 1], [cnt Inf Inf], [stp 1 1]), stafiles, 'uni', 0);
    
    if flag
        htmp = cellfun(@(x) x(:,:,iz,:), htmp, 'uni', 0);
        stmp = cellfun(@(x) x(iz,:,:), stmp, 'uni', 0);
    end
    
    bhis = permute(cat(4, htmp{:}), [1 2 4 3]);
    bsta = permute(cat(3, stmp{:}), [2 3 1]);
    
end
