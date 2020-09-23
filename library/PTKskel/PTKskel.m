% By Ashkan Pakzad on 20th July ashkanpakzad.github.io

% Get skeleton given airway segmentation using PTK without internal loops
% SegImagepath can be any valid image type accepted by PTK e.g. nifti,
% nrrd etc. see github.com/tomdoel/pulmonarytoolkit for details

function PTKskel(SegImagePath, PTK_source)
% SegImagePath is the path to the segmentation image
% PTK_source is the source directory of PTK which need not be given if PTK
% is already in the MATLAB search path.

% add PTK main directory to path if necessary
if nargin > 1
    addpath(PTK_source);
end
PTKAddPaths;
% prevent PTK from creating dialog boxes
ptk_main = PTKMain();

% unzip and load images
if strcmp(SegImagePath(end-2:end), '.gz') == 1
    compressed = 1;
    gunzip(SegImagePath);
    SegImagePath = SegImagePath(1:end-3);
else
    compressed = 0;
end

seg_data = ptk_main.Load(SegImagePath);

skel = seg_data.GetResult('PTKCustomAirwayCentreline');

% save result
[filepath,name,ext] = fileparts(SegImagePath);
savename = [name, '_PTKskel.nii'];

MimSaveNiftiAsNifti(skel, filepath, savename, ptk_main.Reporting);
% compress to .nii.gz
gzip(fullfile(filepath,savename))
delete(fullfile(filepath,savename))

% delete decompression of seg image
if compressed == 1
    delete(fullfile(filepath,[name,ext]))
end
