# Basic use

## Preparing your data
Your data will need to be in Nifti format (.nii or .nii.gz) and read in using `niftiread` and `niftiinfo`. The airway segmentation should be one completely [26-connected](https://en.wikipedia.org/wiki/Pixel_connectivity) object, AQ takes the largest connected object and discards everything else in the segmentation.

It is recommended that the packaged skeletonisation algorithm is used `PTKskel`. The resultant skeleton of this algorithm is suited for use with AirQuant.

## Initialisation
Setting up the AQ object for a case is straightforward. It requires loading up the CT, metadata, airway segmentation and skeleton from [NIFTI format](https://brainder.org/2012/09/23/the-nifti-file-format/). The metadata variable is the output from `niftiinfo`.

Inputs:
* Original CT
* niftiinfo() header information
* Airway Segmentation
* Airway Centreline/Skeleton*
* filename to save results

Output:
* AirQuant object class

*Notes*

Initialising the AQ class calls for the skeleton to be parsed into a directional-graph such that edges face from proximal to the distal direction. [see etc.] This forms the backbone of the AQ class that allows each airway segment to be individually analysed and processed. In addition nodes are bifurcation (or even trifurcation points).

Lobes are classified by examining the graph nodes and their relative position, e.g. the node at the end of the right major bronchi will be more right than the node at the end of the left major bronchi by definition. This is largely based on the method described by [Gu et al.](doi.org/10.1155/2012/382806) with the addition of classifying the Left lingular separately. This can be generated graphically, see [Visualisation](/docs/vis.md).

Each segment's generation is identified by counting the number of nodes from the carina node to origin node of an airway segment. The carina node is first identified by a centrality metric of the first pass to convert the skeleton into a digraph.

All data is converted into LPS orientation (standard for DICOM) once AQ is initialised. So data inside the AQ object may not align with your raw data outside of AQ. It is necessary to declare an orientation within AQ, most CT data is stored in this orientation already. See [Understanding 3D medical image orientation for programmers](https://medium.com/@ashkanpakzad/understanding-3d-medical-image-orientation-for-programmers-fcf79c7beed0).

Functions that take longer to run often call the `save` method, saving the current object's state to a ".mat" file as declared by the `savename` property.

For more information see [library/AirQuant.m](library/AirQuant.m) > methods > %%INITIALISATION > AirQuant.

*Example*
```
% names input files, give fullpaths if not in matlab path
% add AirQuant library to path
AirQuantDir = AirQuantAddPath();
casename = 'github_demo';
dataset = 'example'
results_dir = fullfile(AirQuantDir,'results',dataset, casename);

% Get filenames
CT_name = [casename, '_raw.nii.gz'];
seg_name = [casename, '_seg.nii.gz'];
skel_name = [casename, '_seg_PTKskel.nii.gz'];

% Load CT data as double
meta = niftiinfo(CT_name);
CT = double(niftiread(meta));

% Load Airway segmentation and its skeleton as logicals
S = logical(niftiread(seg_name));
skel = logical(niftiread(skel_name));

% Initialise AirQuant
% Parses CT, segmentation and skeleton to compute airway tree graph
% savename is given to automatically save/load results.
savename = fullfile(results_dir, [casename, '_AQ.mat']);
AQ = AirQuant(CT, meta, S, skel, savename);
```
## Reloading
Previously saved AirQuant objects can be reloaded by calling the AirQuant class with the savename/path only.

*Example*
```
% add AirQuant library to path
AirQuantDir = AirQuantAddPath();
dataset = 'example'
casename = 'github_demo';
results_dir = fullfile(AirQuantDir,'results', dataset, casename);

% Load AirQuant object
savename = fullfile(results_dir, [char(casename, '_AQ.mat']);
AQ = AirQuant(savename);
```

## Batch Processing
AirQuant comes packed with an extended batch function/script under [library/runAQ.m](library/runAQ.m).
This function is designed to run AirQuant on several cases by providing a configuration structure.
It will also run certain visualisation methods and save them to .fig and .png formats.
For alternative batch scripting, it is recommended that you use this function as a template.
An example configuration script can be found under [scripts/example_config.m](scripts/example_config.m).

*Notes*

* This requires all data to be stored exactly within `AirQuant/data/datasetname`.
* The command line output of each case is saved in its own log file under  `AirQuant/data/datasetname/casename/casename_log.txt`.
* If any of the CT, segmentation or skeletonisation datafiles do not exist, it will inform and skip but not throw an error. This is to prevent time wasted in extended jobs.
* If the AQ object already exists in results, it will not skip but instead reload and process it. The AQ traversing method regularly saves the AQ object. If the job is interrupted it just picks up from where it left off.

*Example*
```
AirQuantDir = AirQuantAddPath();

%%% set up the config structure
config = [];
% must be string array, using double dash quotes.
config.casenames = ["github_demo"]; % required
config.dataset = 'example';

% suffix of each file (these are defaults if not provided)
config.CTsuf = '_raw.nii.gz';
config.segsuf= '_seg.nii.gz';
config.skelsuf = '_seg_PTKskel.nii.gz';

%%% pass to AirQuant runner
runAQ(config);

```
