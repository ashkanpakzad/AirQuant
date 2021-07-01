% By Ashkan Pakzad, 2021. ashkanpakzad.github.io
    % batch compute skeleton of input airway graph using PTKskeleton.
    % This method ensures there are no internal loops within the skeleton.
    
% This batch script assumes that all of your cases are saved within their
% own folders and have a folder within them that contains the segmentation
% nifti file in compressed nifti format i.e. case_seg.nii.gz.
    % parentdatapath
    %       |- case
    %           |- data_folder
    %               |- case_seg.nii.gz
    %               |- case_seg_PTKskel.nii.gz (generated)  
% The result along with log file and preview visualisation will be saved in 
% the same directory as the segmentation.
% 
AirQuantAddPath(); % Ensure all AirQuant files are in matlab path

%% Set up
tic
% path to data folder
parentdatapath = 'C:\data'; 
data_folder = 'postprocessed'; % folder within casename folder that contains the segmentation
seg_suffix = '_seg.nii.gz';
skel_suffix = '_PTKskel.nii.gz';

%% run
% Get list of casenames
casedir = dir(parentdatapath);

% function for full path of listing
fullpath = @(x, i) fullfile(x(i).folder, x(i).name);

for ii = 3:length(casedir)
    disp(ii-2)
    casename = casedir(ii).name;
    casepath = fullpath(casedir,ii);
    disp(casename)

    % find segmentation path
    preddir = fullfile(casepath,data_folder); % prediction TIFF files
    savepath = preddir;
    overwritelog = false; % If log should be reset if already exists (false to append)
    
    % check prediction path exists
    if ~isfolder(preddir)
        warning([preddir, 'folder containing segmentation does not exist, SKIPPING.']);
        continue
    end
    
    % set log
    logname = [casename, '_PTKskel_log.txt'];
    diary(fullfile(preddir, logname))
    disp(['Running PTKskel on ', casename]);
    disp(datetime)
    
    % change directory (current folder) to case prediction
    cd(preddir)
    
    % check segmentation exists
    segname = [casename, seg_suffix];
    if ~isfile(segname)
        warning ([segname, ' does not exist, SKIPPING (ensure the segmentation exists to run PTKskel).'])
        continue
    end
    
    % check if skeleton already exists
    % assumes gzip compression
    skelname = [segname(1:end-7), skel_suffix];
    if isfile(skelname)
        warning ([skelname, ', skeleton already exists, SKIPPING (delete to run PTK again).'])
        continue
    end
    
    % RUN PTKSKEL
    % must be in matlab current path
    PTKskel(segname, 1);

    % reload to ensure read image is suitable
    seg = niftiread(segname);
    skel = niftiread(skelname);
    
    % visualise and save
    fig = figure;
    patch(isosurface(seg),'EdgeColor', 'none','FaceAlpha',0.3);
    title(casename)
    hold on
    patch(isosurface(skel),'EdgeColor', 'r',...
        'FaceAlpha',1);
    isosurface(skel)
    axis vis3d
    view(80,0)
    saveas(fig, [casename,'_PTKskel.png']);
    
    % reset
    close all;
    diary off;
end
disp(toc/60)
