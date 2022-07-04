% By Ashkan Pakzad on 20th July ashkanpakzad.github.io

% Get skeleton given airway segmentation using PTK removing internal loops
% SegImagepath can be any valid image type accepted by PTK e.g. nifti,
% nrrd etc. see github.com/tomdoel/pulmonarytoolkit for details
%
% Seed is 3x1, [x,y,z] vector of coordinates in LPS coordinates. It is
% converted into LPI for compatibility with PTK standard image format. and
% first two axes swapped to image space [row, col].
%

function skelvol = PTKskel(segpath, seed, saveflag)
% SegImagePath is the path to the segmentation image
% PTK_source is the source directory of PTK which need not be given if PTK
% is already in the MATLAB search path.
% add PTK main directory to path if necessary

if nargin < 2
    seed = [];
end

if nargin < 3
    saveflag = 1;
end

% init PTK
ptk_main = PTKMain();

% get full segpath and load
segpath_full = which(segpath);
seg_data = ptk_main.Load(segpath_full);
seg_ptk = seg_data.GetResult('PTKOriginalImage');

if isempty(seed) 
% get seed by trachea identification if not specified
    [seed, ~] = PTKFindTopOfTrachea(seg_ptk, ptk_main.Reporting, false);
else
    % swap axis 1 and 2 and convert to LPI
    imsz = seg_ptk.OriginalImageSize;
    seed = [seed(2), seed(1), imsz(3) - seed(3) + 1];
end

skel_ptk = SeededCentreline(ptk_main, seg_ptk, seed);

% save result if requested and segmentation loaded originally.
[filepath,name,~] = fileparts(segpath_full);

if endsWith(name, '.nii')
    name = name(1:end-4);
end

if saveflag == 1
    savename = [name, '_PTKskel.nii'];
    MimSaveNiftiAsNifti(skel_ptk, filepath, savename, ptk_main.Reporting);
    % compress to .nii.gz
    gzip(fullfile(filepath,savename))
    delete(fullfile(filepath,savename))    
end

if nargout > 0
    skelvol = skel_ptk.RawImage;
    % from PTK orientation to original
    skelvol = permute(flip(skelvol,3),[2,1,3]);
end
