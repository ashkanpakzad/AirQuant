% By Ashkan Pakzad on 20th July ashkanpakzad.github.io

% Get skeleton given airway segmentation using PTK removing internal loops
% SegImagepath can be any valid image type accepted by PTK e.g. nifti,
% nrrd etc. see github.com/tomdoel/pulmonarytoolkit for details

function skelvol = PTKskel(SegImage, saveflag)
% SegImagePath is the path to the segmentation image
% PTK_source is the source directory of PTK which need not be given if PTK
% is already in the MATLAB search path.
% add PTK main directory to path if necessary
if nargin < 2
    saveflag = 0;
end
% prevent PTK from creating dialog boxes
ptk_main = PTKMain();

% load in segmentation and execute

% load in using niftiread and convert to PTK image
if strcmp(SegImage(end-2:end), '.gz') == 1
    compressed = 1;
    gunzip(SegImage);
    SegImage = SegImage(1:end-3);
else
    compressed = 0;
end

seg_data = ptk_main.Load(SegImage);

skel = seg_data.GetResult('PTKCustomAirwayCentreline');

% save result if requested and segmentation loaded originally.
[filepath,name,ext] = fileparts(SegImage);

if saveflag == 1
    savename = [name, '_PTKskel.nii'];
    
    MimSaveNiftiAsNifti(skel, filepath, savename, ptk_main.Reporting);
    % compress to .nii.gz
    gzip(fullfile(filepath,savename))
    delete(fullfile(filepath,savename))
    
    % delete decompression of seg image
    
end

if compressed == 1
    delete(fullfile(filepath,[name,ext]))
end

if nargout > 0
    skelvol = skel.RawImage;
    % Test output form.
    skelvol = permute(flip(skelvol,3),[2,1,3]);
end
