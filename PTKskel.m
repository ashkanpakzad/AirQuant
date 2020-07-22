% By Ashkan Pakzad on 20th July ashkanpakzad.github.io

% Get skeleton given airway segmentation using PTK without internal loops
% SegImagepath can be any valid image type accepted by PTK e.g. nifti,
% nrrd etc. see github.com/tomdoel/pulmonarytoolkit for details

function result = PTKskel(SegImagePath, PTK_source)
% SegImagePath is the path to the segmentation image
% PTK_source is the source directory of PTK which need not be given if PTK
% is already in the MATLAB search path.

% add PTK main directory to path if necessary
if nargin > 1
    addpath(PTK_source)
end
PTKAddPaths;
% prevent PTK from creating dialog boxes
reporting = CoreReporting([], false);
ptk_main = PTKMain(reporting);

% unzip and load images
if strcmp(SegImagePath(end-3:end), '.gz') == 1
    gunzip(SegImagePath)
    SegImagePath = SegImagePath(1:end-3);
end

seg_data = ptk_main.Load(SegImagePath);

% get image
skel = seg_data.GetResult('PTKCustomAirwayCentreline');
result = skel.RawImage;
