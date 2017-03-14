function Data = ncreads(file, varargin)
%NCVARSGET Read several variables from a netcdf file
%
% Data = ncreads(file, var1, var2, ...)
% Data = ncreads(file, Scs, var1, var2, ...)
%
% This function reads multiple variables from a netcdf file.  It is just a
% wrapper function for ncread.
%
% Input variables:
%
%   file:   name of netcdf file
%
%   var:    name of variable to read in, which must correspond exactly to a
%           variable in the file.  If no variable names are listed, all
%           variables are read.  Can also include the string 'dimensions'
%           to return variables whose name and size match one of the
%           dimensions of the file.
%
%   Scs:    1 x 1 structure of start-count-stride values.  Field names
%           should match names of dimension variables in the file.  If
%           provided, any variable with that dimension will be subset as
%           described.
%         
% Output variables:
%
%   Data:   1 x 1 structure, with fieldnames corresponding to requested
%           variables, each holding the data value for that variable

% Copyright 2015-2017 Kelly Kearney

% Parse input

if ~exist(file, 'file')
    error('File not found');
end

isscs = cellfun(@isstruct, varargin);
if any(isscs)
    Scs = varargin{isscs};
    varnames = varargin(~isscs);
else
    Scs = struct;
    varnames = varargin;
end

Info = ncinfo(file);
if isempty(varnames)
    varnames = {Info.Variables.Name};
end

scsfld = fieldnames(Scs);
[tf,loc] = ismember(scsfld, {Info.Dimensions.Name});
if ~all(tf)
    str = sprintf('%s,', scsfld{~tf});
    error('Field in Scs input (%s) did not match file dimension names', str(1:end-1));
end

% Set up start-stride-count for all dimensions

nd = length(Info.Dimensions);

start = ones(1,nd);
count = ones(1,nd)*Inf;
stride = ones(1,nd);

if ~isempty(loc)
    scs = struct2cell(Scs);
    scs = cat(1, scs{:});
    start(loc)  = scs(:,1);
    count(loc)  = scs(:,2);
    stride(loc) = scs(:,3);
end

% If 'dimensions' is passed as variable name, look for variables that share
% name and size with a dimension

isdimshortcut = strcmp(varnames, 'dimensions');
if any(isdimshortcut)
    varnames = varnames(~isdimshortcut);
    if ~exist('Info', 'var')
        Info = ncinfo(file);
    end
    filevars = {Info.Variables.Name};
    if any(strcmp('dimensions', filevars))
        warning('''dimensions'' refers to a variable in this file; option disabled');
    else
        dims = {Info.Dimensions.Name};
        isvar = ismember(dims, filevars);
        dimvars = dims(isvar);
        varnames = [dimvars varnames];
        varnames = unique(varnames);
    end
end

% Locate and check variable names

nvar = length(varnames);

[vfound, vloc] = ismember(varnames, {Info.Variables.Name});
if ~all(vfound)
    str = sprintf('%s,', varnames{~vfound});
    error('Variable names (%s) not found in file', str(1:end-1));
end

% Read data

for iv = 1:nvar
    if isempty(Info.Variables(vloc(iv)).Dimensions)
        Data.(varnames{iv}) = ncread(file, varnames{iv});
    else
        [~,dloc] = ismember({Info.Variables(vloc(iv)).Dimensions.Name}, {Info.Dimensions.Name});
        Data.(varnames{iv}) = ncread(file, varnames{iv}, start(dloc), count(dloc), stride(dloc));
    end
end