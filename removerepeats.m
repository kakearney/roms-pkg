function [A, x] = removerepeats(A, var, files, varargin)

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

