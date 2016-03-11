function Data = ncreads(file, varargin)
%NCVARSGET Read several variables from a netcdf file
%
% Data = ncvarsget(file, var1, var2, ...)
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
%           variables are read.
%         
% Output variables:
%
%   Data:   1 x 1 structure, with fieldnames corresponding to requested
%           variables, each holding the data value for that variable

% Copyright 2015 Kelly Kearney

varnames = varargin;
if isempty(varnames)
    Info = ncinfo(file);
    varnames = {Info.Variables.Name};
end

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

nvar = length(varnames);

for iv = 1:nvar
    Data.(varnames{iv}) = ncread(file, varnames{iv});
end