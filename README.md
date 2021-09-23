
# Some ROMS utilities


Author: Kelly Kearney


This repository holds a loose collections of utilities I've written to interact with output from the [Regional Ocean Modeling System (ROMS)](https://www.myroms.org), a widely-used primitive-equations ocean model.



## Contents

            
- Getting started        
- Functions        
- Contributions

## Getting started


**Prerequisites**


For full functionality, these tools require Matlab R2012a or later with the Mapping, Statistics, and Image Processing Toolboxes.


**Downloading and installation**


This code can be downloaded from [Github](https://github.com/kakearney/roms-pkg/).


**Matlab Search Path**


The following folders need to be added to your Matlab Search path (via `addpath`, `pathtool`, etc.):



```
roms-pkg/roms
roms-pkg/cellstr2
roms-pkg/minmax
roms-pkg/regexpfound
roms-pkg/setgetpos_V1.2
```


The primary functions for this toolbox can be found in the roms folder; the remaining folders hold dependent functions used by thse primary ones.


In addition to the dependencies listed above, some roms-pkg functions rely on the netCDF-reading utilities from the [Climate Data Toolbox](https://github.com/chadagreene/CDT); I highly recommend users download this entire toolbox for its extra netCDF utilities and its variety of plotting and analysis functions that can aid in working with ROMS output datasets.  Finally, the ROMS folder includes two utility functions, `stretching` and `set_depth`, that are copied over from the [myroms.org Matlab Scripts](https://www.myroms.org/wiki/Matlab_Scripts); personally, I'm not a fan of this toolbox due to its clunky syntax, but it's useful to have handy occasionally if you work with ROMS output.



## Functions



  - `plotromsrho`: Plot a 2D slice of a rho variable to projected pcolor   plot
  - `animateromsrho`: Animate a 2D slice of a rho variable over time
  - `calcromsz`: Calculate rho- and w-depth values based on bottom depth   and surface height
  - `romsboundarycoords`: Extract boundary slices from 4D grid
  - `romsgeometryparams`: Calculate various grid-related geometry   parameters
  - `parsevarinfo`: Read data from a varinfo.dat file into a table
  - `parseromslog`: Read archived standard output text from from a ROMS   simulation


## Contributions


Community contributions to this package are welcome!


To report bugs, please submit [an issue](https://github.com/kakearney/example-pkg/issues) on GitHub and include:



  - your operating system
  - your version of Matlab and all relevant toolboxes (type `ver` at the Matlab command line to get this info)
  - code/data to reproduce the error or buggy behavior, and the full text of any error messages received

Please also feel free to submit enhancement requests, or to send pull requests (via GitHub) for bug fixes or new features.


I do monitor the MatlabCentral FileExchange entry for any issues raised in the comments, but would prefer to track issues on GitHub.



<sub>[Published with MATLAB R2019a]("http://www.mathworks.com/products/matlab/")</sub>