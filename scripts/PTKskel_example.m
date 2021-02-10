% By Ashkan Pakzad, 2020. ashkanpakzad.github.io
    % Compute skeleton of input airway graph using PTKskeleton.
    % This method ensures there are no internal loops within the skeleton.
    % Dependent on PTK github
    
segname = 'github_demo_seg.nii.gz';
AirQuantAddPath(); % Ensure all AirQuant files are in matlab path

% must be in matlab current path
PTKskel(segname, 1);

%% check
% assumes gzip compression
skelname = [segname(1:end-7), '_PTKskel.nii.gz'];

seg = niftiread(segname);
skel = niftiread(skelname);

figure;
patch(isosurface(seg),'EdgeColor', 'none','FaceAlpha',0.3);
hold on
isosurface(skel)
vol3daxes(obj)