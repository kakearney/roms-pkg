function [xpoly,ypoly] = mask2poly(mx, my, mask)
%MASK2POLY Calculates polygon outline of a gridded mask
%
% [x,y] = mask2poly(mx, my, mask)
%
% This function finds the outline of the polygons defined by a gridded
% mask, similar to bwboundaries except that the grid is defined by its edge
% coordinates in pcolor-style. This requires that the edge definitions are
% one cell longer than the mask itself.
%
% Input variables:
%
%   mx:     n x 1 or m x n array defining the edges of each cell in the x
%           direction
%
%   my:     m x 1 or m x n array defining the edges of each cell in the y
%           direction
%
%   mask:   m-1 x n-1 array, with 1s indicating the mask region and 0 the
%           background
%
% Output variables:
%
%   xpoly:  x coordinates of mask polygon, NaN-delimted between segments
%
%   ypoly:  y coordinates of mask polygon, NaN-delimted between segments

% Copyright 2013 Kelly Kearney

% Check input

if isvector(mx) & isvector(my)
    [mx, my] = meshgrid(mx, my);
end

[nrow, ncol] = size(mask);

if ~isequal(size(mx), size(my), [nrow+1 ncol+1])
    error('mx and my must be one row and column larger than mask');
end

% Edges of each cell

[xv, yv] = deal(zeros(nrow,ncol,5));
for ir = 1:nrow
    for ic = 1:ncol
        xv(ir,ic,:) = [mx(ir,ic) mx(ir,ic+1) mx(ir+1,ic+1) mx(ir+1,ic) mx(ir,ic)]';
        yv(ir,ic,:) = [my(ir,ic) my(ir,ic+1) my(ir+1,ic+1) my(ir+1,ic) my(ir,ic)]';
    end
end

xv = reshape(xv, [], 5);
yv = reshape(yv, [], 5);

idx = find(mask);
xedge = xv(idx,:);
yedge = yv(idx,:);


xyedge = [...
    xedge(:,1:2) yedge(:,1:2)
    xedge(:,2:3) yedge(:,2:3)
    xedge(:,3:4) yedge(:,3:4)
    xedge(:,4:5) yedge(:,4:5)];

% Orient all edges the same way, so tracing direction doesn't
% matter

isokay = xyedge(:,1) < xyedge(:,2) | ...
        (xyedge(:,1) == xyedge(:,2) & xyedge(:,3) < xyedge(:,4));

xyedge(~isokay,:) = xyedge(~isokay, [2 1 4 3]);

% Locate the edges only repeated once; these are the perimeter edges

[unqedge, nedge] = consolidator(xyedge, [], 'count');
xunq = unqedge(nedge==1, 1:2);
yunq = unqedge(nedge==1, 3:4);

xunq = [xunq'; nan(1,size(xunq,1))];
yunq = [yunq'; nan(1,size(yunq,1))];

% Merge line segments into one or more polygons

[xpoly, ypoly]  = polymerge(xunq(:), yunq(:));






