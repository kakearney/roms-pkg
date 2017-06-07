function [A, x] = removerepeats(A, var, files, varargin)
%REMOVEREPEATS Removes non-unique values along one dimensions
%
% [B, x] = removerepeats(A, var, files, ...)
%
% Ths function is intended as a companion function to ncreadsseries.  It
% allows one to remove any repeated values along a dimension, as determined
% by the values in a one-dimensional variable along that dimension.  
%
% My original intent was to add this functionality to ncreadsseries itself,
% but this is difficult due to the flexibility of netCDF files, where
% dimensions don't necessarily include a one-to-one match with a particular
% variable.  This function is primarily intended to work along the
% unlimited dimension along which ncreadsseries data was concatenated, and
% has not been throroughly.  
%
% Input variables:
%
%   A:      structure of netCDF data (output of ncreads or ncreadsseries,
%           or similar) 
%
%   var:    name of variable a one-dimensional variable in the file on
%           which the unique calculation will be based.  Data in A along
%           this dimension should not be subset (i.e. should not have been 
%           read in with a start-count-stride modifier).
%
%   files:  files from which the data in A were read.
%
% Optional input variables (passed as parameter value pairs):
%
%   pos:    position (occurrence) to keep in data ('first' or 'last').  
%           Default: 'last'
%
% Output variables:
%
%   B:      structure of netCDF data, matching A but modified to remove
%           data corresponding to repeats along the indicated dimension.
%           Variables that do not include that dimension will be
%           unmodified.
%
%   x:      values of var, with repeats removed.

% Copyright 2017 Kelly Kearney

p = inputParser;
p.addParameter('pos', 'last', @(x) validateattributes(x, {'char','string'}));
p.parse(varargin{:});

Opt = p.Results;


if ischar(files)
    Tmp = ncreads(files, var);
    reffile = files;
elseif iscellstr(files)
    Tmp = ncreadsseries(files, var);
    reffile = files{1};
end

% Dimension corresponding to variable

I = ncinfo(reffile, var);
if length(I.Dimensions) > 1
    error('Reference variable must be one-dimensional');
end
dimname = I.Dimensions.Name;

% Find unique indices of reference variable

[x, ix] = unique(Tmp.(var), Opt.pos);

flds = fieldnames(A);
for iv = 1:length(flds)
    Vinfo = ncinfo(reffile, flds{iv});
    if ~isempty(Vinfo.Dimensions)
        [tf, loc] = ismember(dimname, {Vinfo.Dimensions.Name});
        if tf
            nd = length(Vinfo.Dimensions);
            if nd == 1
                A.(flds{iv}) = A.(flds{iv})(ix);
            else
                neworder = [loc setdiff(1:nd, loc, 'stable')];
    
                tmp = permute(A.(flds{iv}), neworder);
                sz = size(tmp);
                tmp = reshape(tmp, sz(1), []);
                tmp = tmp(ix,:);
                tmp = reshape(tmp, [length(x), sz(2:end)]);

                A.(flds{iv}) = ipermute(tmp, neworder);
            end
        end
    end   
end

