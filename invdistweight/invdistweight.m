function [u, dists] = invdistweight(varargin)
%INVDISTWEIGHT Inverse distance weighting interpolation
%
% u = invdistweight(uobs, p, xy, xyobs)
% u = invdistweight(uobs, p, xy, xyobs, coord)
% u = invdistweight(uobs, p, 'dist', dists)
% [u, dist] = invdistweight(...)
% 
% This function performs inverse distance weighting interpolation based on
% known values for a set of scattered points, using either Euclidean or
% geographic distance metrics.
%
% Input variables:
%
%   uobs:   n1 x 1 vector, values at observed points
%
%   p:      scalar, >0, power parameter for weight calculations 
%           w_ij = 1./(d_ij^p), where d_ij is the distance between points i
%           and j
%   
%   xy:     n2 x 2 matrix, x and y coordinates where interpolated values
%           will be calculated
%   
%   xyobs:  n1 x 2 matrix, x and y coordinates of observed points
%
%   coord:  'cart' or 'geo', specifying whether coordinates are defined in
%           cartesion (x/y) or geographic (lon/lat) coordinates.
%
%   dists:  n1 x n2 distance matrix.  Include in place of xy and xyobs in
%           order to provide your own distance metric, where dists(i,j) is
%           the distance from observation i to interpolation point j.  In
%           addition to being useful for flexible definition of distance,
%           this can also speed up calculations when using the same
%           observation and destination points with different datasets,
%           since the distance calculation is the most time-consuming part
%           of this code. 
%
% Output variables:
%
%   u:      n2 x 1 vector, interpolated values
%
%   dist:   n1 x n2 distance matrix, where dist(i,j) is the distance from
%           observation i to interpolation point j.  For geographic
%           coordinates, units are km.

% Copyright 2015 Kelly Kearney

%-----------------------
% Parse input
%-----------------------

% Parsers

p1 = inputParser;
p1.addRequired('uobs',   @(x) validateattributes(x, {'numeric'}, {'column'}));
p1.addRequired('p',      @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p1.addRequired('xy',     @(x) validateattributes(x, {'numeric'}, {'ncols', 2}));
p1.addRequired('xyobs',  @(x) validateattributes(x, {'numeric'}, {'ncols', 2}));
% p1.addOptional('xypoly', [],     @(x) validateattributes(x, {'numeric'}, {'ncols', 2}));
p1.addOptional('coord',  'cart', @(x) validateattributes(x, {'char'}, {}));

p2 = inputParser;
p2.addRequired('uobs',   @(x) validateattributes(x, {'numeric'}, {'column'}));
p2.addRequired('p',      @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p2.addParameter('dist', [],      @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Make sure input fits one of the two options

try
    parse(p1, varargin{:})
    In = p1.Results;
catch ME
    parse(p2, varargin{:});
    In = p2.Results;
end

% Check and expand string if necessary

if isfield(In, 'coord')
    In.coord = validatestring(In.coord, {'cart', 'geo'}, '', 'coord');
end

%-----------------------
% Interpolate
%-----------------------

% Calculate distances between observations and interpolation points

if isfield(In, 'dist')    % No calculation necessary
    dists = In.dist;
    if size(dists,1) ~= length(In.uobs)
        error('dist must have same number of rows as uobs');
    end
else
    
    nobs = size(In.xyobs,1);
    npt  = size(In.xy,1);
    
    [srcidx, snkidx] = ndgrid(1:nobs, 1:npt);
    xy1 = In.xyobs(srcidx,:);
    xy2 = In.xy(snkidx,:);
    switch In.coord
        case 'cart'
            dists = sqrt(sum((xy1 - xy2).^2,2));
            dists = reshape(dists, size(srcidx));
        case 'geo'
            dists = deg2km(distance(xy1(:,2), xy1(:,1), xy2(:,2), xy2(:,1))); 
            dists = reshape(dists, size(srcidx));
    end

end
    
% Inverse distance weighting

if any(isnan(In.uobs))
    warning('IDW:NaN', 'Ignoring NaNs');
    isn = isnan(In.uobs);
    In.uobs = In.uobs(~isn);
    dists = dists(~isn,:);
end

        
w = 1./(dists.^In.p);
wu = bsxfun(@times, w, In.uobs);

wsum = sum(w, 1);

u = sum(bsxfun(@rdivide, wu, wsum), 1)';










