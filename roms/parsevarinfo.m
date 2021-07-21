function vtbl = parsevarinfo(vfile)
%PARSEVARINFO Read varinfo.dat data into a table
%
% vtbl = parsevarinfo(vfile)
%
% This function reads in the elements of a ROMS varinfo.dat file and
% formats them as a table.
%
% Input variables:
%
%   vfile:  path to varinfo.dat file
%
% Output variables:
%
%   vtbl:   table with the following variables, all variables are character
%           arrays:
%
%           fieldname:  field variable name, Vinfo(1)
%
%           longname:   long-name attribute, Vinfo(2)
%
%           units:      units attribute, Vinfo(3)
%
%           type:       field type attribute, Vinfo(4)
%
%           timevar:    associated time variabel name, Vinfo(5)
%
%           index:      index variable name used in information arrays,
%                       Vinfo(6)
%
%           cgrid:      staggered C-grid variable type, Vinfo(7)
%
%           scale:      scale to convert input data to model units,
%                       Vinfo(8)
%
%           x_comment:  any commented text found to the right of each of
%                       the above fields

% Copyright 2021 Kelly Kearney

vtext = fileread(vfile);
vtext = regexp(vtext, '\n', 'split');

% Remove comment lines, blank lines, and URL line
allcomment = regexpfound(vtext, '^\s*!');
isemp = cellfun(@isempty, vtext) | regexpfound(vtext, '^\s*$');
isurl = startsWith(vtext, '$URL') | startsWith(vtext, '''$URL');

vtext = vtext(~isemp & ~allcomment & ~isurl);

vtext = reshape(vtext, 8,[]);
vtext = regexprep(vtext, '^\s*''', '');
vtext = regexprep(vtext, '''\s*$', '');

vtbl = makevtable(vtext');

% Parse out comments

fld = vtbl.Properties.VariableNames;
for ii = 1:length(fld)
    cfld = [fld{ii} '_comment'];
    Comment.(cfld) = unpackcell(regexp(vtbl.(fld{ii}), '!.*', 'match'));
    Comment.(cfld) = regexprep(Comment.(cfld), '\s*!\s*', '');

    vtbl.(fld{ii}) = strtrim(regexprep(vtbl.(fld{ii}), '!.*', ''));
    vtbl.(fld{ii}) = regexprep(vtbl.(fld{ii}), '''\s*$', '');
end

vtbl = [vtbl struct2table(Comment)];
    
end

% Subfunctions

function v = makevtable(x)
    v = cell2table(x, 'VariableNames', ...
        {'fieldname', 'longname', 'units', 'type', 'timevar', 'index', 'cgrid', 'scale'});
end

function x = unpackcell(x)

    for ii = 1:length(x)
        if isempty(x{ii})
            x{ii} = '';
        else
            x{ii} = x{ii}{1};
        end
    end
end