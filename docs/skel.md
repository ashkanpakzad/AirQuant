# Skeletonisation

The airway skeleton/centreline should be true to the natural centreline of the airway and contain no loops such that it can be broken down into a tree object. AQ parses the skeleton into a network digraph allowing analysis on individual airway segments.

AQ will warn of any anomalies in the skeleton on initialisation and flag them on visualisation of the skeleton graph but make no effort to solve them.

A robust skeleton is necessary for successful analysis with AirQuant, to this end a suitable algorithm, [library/PTKskel] is packaged with AirQuant. Example use can also be found in [scripts/example_skel]

## PTKskel
Provide the filename of the segmentation to be skeletonised. The result is saved in the same folder and the same name with an appended "_PTKskel". The saved nifti file will be in alignment with the original segmentation file, allowing it to be used easily outside of PTK/AirQuant.

The PTK library requires that the input segmentation be in the MATLAB current folder.

*Notes*

This algorithm utilises the library of the [PulmonaryToolKit (PTK) by Tom Doel](https://github.com/tomdoel/pulmonarytoolkit) to skeletonise an already complete segmentation. Original PTK plugins only allows skeletonisation of airway segmentations complete by its own algorithm, this algorithm essentially calls PTK region growing algorithm to re-segment the airways using a propagating wavefront method originating from the trachea, this algorithm parses the airways at the same time. The resultant PTKairways object can then be passed to the PTKskeletonisation library to employ the algorithm based on [Palágyi et al.](doi.org/10.1016/j.compbiomed.2005.05.004). It also checks for loops and removes 'offending' skeleton branches.

AirQuant has methods to plot the segment and skeleton in one figure, see [Visualisation](/docs/vis.md).

*Example*
```
% must be in matlab current path
segname = 'github_demo_seg.nii.gz';
% Ensure all AirQuant files are in matlab path
AirQuantAddPath();

PTKskel(segname);
```
