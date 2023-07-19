function yi = fastinterpcol(x,y,xi)
%INTERPCOL Quick 1-D linear interpolation along the column dimension.
%
%   YI = FASTINTERPCOL(X,Y,XI) returns linear interpolation result along the
%   column dimension.
%   X, Y, XI must have the same number of rows. Every column of X must be
%   monotonically increasing.
%   This source code is an expanded version of interp1q. It uses
%   vectorization to speed up the calculation over the column dimension.
%
%   Modified from INTERP1q to allow vectorization of the linear
%   interpolation along the column dimension.
%   Tan H. Nguyen, MIT, 2019.
%   For question: please send to thnguyn@mit.edu.

    % Step 1: find the indices of the xi's elements in a combined [x; xc] sorted array.
    xc = [x; xi];
    [~, xcIdx] = sort(xc, 1);  % Find xcIdx such that sortedXc = xc(xcIdx);
    colIdx1d = 1 : size(x, 2); % Generate a column index vector
    colIdx = repmat(colIdx1d, [size(xc, 1) 1]);
    xcIdxCombined = sub2ind(size(xc), xcIdx, colIdx);  %Convert to a 1D coordinate vector
    
    r = zeros(size(xc));
    linearIdx = 1 : size(xc, 1);
    linearIdx = repmat(linearIdx', [1 size(xc, 2)]);
    r(xcIdxCombined) = linearIdx;         % Find r such that for every column, xc = sortedXc(r);
    rXiIdx = r(size(x, 1) + 1 : end, :);  % Extract the index of the xi in the sorted array.

    % Step2: Find the index of the nearest neighbor in x that is closest to every
    % each element in Xi. This is the index in the sorted x array,
    offsetIdx = 1 : size(rXiIdx, 1);
    lhsXIdx = rXiIdx - repmat(offsetIdx', [1 size(rXiIdx, 2)]);
    % Make sure that we can handle the right boundary well by never have the
    % leftIdx pointed to the last element.
    lhsXIdx(lhsXIdx == size(x, 1)) = size(x, 1) - 1;

    % Step 3: interpolation.
    colIdx = repmat(colIdx1d, [size(xi, 1) 1]);
    lhsXIdxCombined = sub2ind(size(x), lhsXIdx, colIdx);
    xlhs = x(lhsXIdxCombined);
    ylhs = y(lhsXIdxCombined);
   
    rhsXIdxCombined = sub2ind(size(x), lhsXIdx + 1, colIdx);
    xrhs = x(rhsXIdxCombined);
    yrhs = y(rhsXIdxCombined);
    u = (xi - xlhs) ./ (xrhs - xlhs);
    yi = (1 - u) .* ylhs + u .* yrhs;
end
