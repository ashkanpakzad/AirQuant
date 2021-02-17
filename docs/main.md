# AirQuant Docs
By Ashkan Pakzad https://ashkanpakzad.github.io

See [/readme.md](/readme.md) before reading this document regarding install, about and purpose.


This software primarily revolves around the [library/AirQuant.m](library/AirQuant.m) class. Every CT case will exist as a new AirQuant class object. When a case is loaded in, it is first initialised, performing a number of checks and short processes.
Currently there is only one method for taking airway measurements, the [FWHMesl method](https://doi.org/10.1117/12.595283). This first relies on the CT scan to be interpolated at right angles to the principle airway axis. For scientific reading on the concept, please see  [K Quan et al.](https://doi.org/10.1117/12.2292306) on which this software is based.

As most of the code behind this software is in [library/AirQuant.m](library/AirQuant.m) it is recommended that you explore it in MATLAB by collapsing all code and expanding as necessary. Set code collapsing for section

Example use can also be found in [scripts/example.m](scripts/example.m) and example data is in [data/](data/).




# Contents:
 * [Data and Initialisation](/docs/basic.md) - prepping data and initialisation an AirQuant object.
 * [PTKskel: Skeletonisation](/docs/skel.md) - Skeletonisation and the packaged algorithm.
 * [CT Airway Interpolation](/docs/interp.md) - Interpolation of the CT along the airway principal axis.
 * [FWHMesl method](/docs/fwhm.md) - Measuring the airway.
 * [Tapering Metrics](/docs/taper.md) - Tapering metrics that can be generated.
 * [Visualisation](/docs/vis.md) - Methods for visualising global data.
 * [Segmental Visualisation](/docs/segvis.md) - Methods for visualising data for individual segments.

 * [Future Ambitions](/docs/future.md) - watch this space.
