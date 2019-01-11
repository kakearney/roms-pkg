function u = idw(uobs, d, p, dmax)
%IDW Inverse distance weighting calc
%
% This function performs an inverse distance weighting interpolation,
% assuming all distance values have been pre-caculated.  
%
% u = idw(uobs, d, p)
% u = idw(uobs, d, p, dmax)
%
% Input variables:
%
%   uobs:   n x m ... x npt array, values at closest observation points to
%           the n x m x ... array being estimated.
%
%   d:      n x m ... x npt array, distance to closest observation points
%
%   p:      power parameter
%
% Output variables:
%
%   u:      n x m x ... array, interpolated value

% Copyright 2017 Kelly Kearney

nd = ndims(uobs);

if nargin > 3
    w = (max(0, dmax - d)./(dmax .* d)).^(1/p);
else
    w = 1./d.^p;
end
    
wsum = nansum(w,nd);

wu = w .* uobs;
u = nansum(wu,nd)./wsum;

isn = all(isnan(uobs),nd);
u(isn) = NaN;


