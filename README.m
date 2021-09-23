%% Some ROMS utilities
% Author: Kelly Kearney
%
% This repository holds a loose collections of utilities I've written to
% interact with output from the <https://www.myroms.org Regional Ocean
% Modeling System (ROMS)>, a widely-used primitive-equations ocean model.
% The primary functions can be found in the
%
% Paragraph description for this function or suite of functions.
%
%% Getting started
%
% *Prerequisites*
%
% This function requires Matlab R14 or later.
%
% *Downloading and installation*
%
% This code can be downloaded from <https://github.com/kakearney/example-pkg/ Github>
% or the
% <http://www.mathworks.com/matlabcentral/fileexchange/xxxx-example
% MatlabCentral File Exchange>.  The File Exchange entry is updated daily
% from the GitHub repository.   
%
% *Matlab Search Path*
%
% The following folders need to be added to your Matlab Search path (via
% |addpath|, |pathtool|, etc.):
%
%   roms-pkg/roms
%   roms-pkg/cellstr2
%   roms-pkg/minmax
%   roms-pkg/regexpfound
%   roms-pkg/setgetpos_V1.2
%
% The primary functions for this toolbox can be found in the roms folder,
% with additional child functions called by these within the remaining
% folders.

%% Functions
%
% * 
%
%  y = example(x)
%
% Input variables:
%
% * 
% * 

%% Examples

% Here's some example code, with an image

x = 1:10;
y = x.^2

plot(x,y);

%% Contributions
%
% Community contributions to this package are welcome!
% 
% To report bugs, please submit
% <https://github.com/kakearney/example-pkg/issues an issue> on GitHub and
% include:  
% 
% * your operating system
% * your version of Matlab and all relevant toolboxes (type |ver| at the Matlab command line to get this info)  
% * code/data to reproduce the error or buggy behavior, and the full text of any error messages received 
% 
% Please also feel free to submit enhancement requests, or to send pull
% requests (via GitHub) for bug fixes or new features. 
% 
% I do monitor the MatlabCentral FileExchange entry for any issues raised
% in the comments, but would prefer to track issues on GitHub. 
% 

