% By Ashkan Pakzad, 2020. ashkanpakzad.github.io

%% Set-up directories
% names input files, give fullpaths if not in matlab path
% add AirQuant library to path
AirQuantDir = AirQuantAddPath();
casename = 'github_demo';
results_dir = fullfile(AirQuantDir,'results', casename);

%% Get filenames
CT_name = [casename, '_raw.nii.gz'];
seg_name = [casename, '_seg.nii.gz'];
skel_name = [casename, '_seg_PTKskel.nii.gz'];

% Load CT data as double
meta = niftiinfo(CT_name);
CT = double(niftiread(meta));

% Load Airway segmentation and its skeleton as logicals
S = logical(niftiread(seg_name));
skel = logical(niftiread(skel_name));

%% Initialise AirQuant
% Parses CT, segmentation and skeleton to compute airway tree graph
% savename is given to automatically save/load results.
savename = fullfile(results_dir, [casename, '_AQ.mat']);
AQ = AirQuant(CT, meta, S, skel, savename);

%% Display Airway Graph in 3D and 2D
figure
PlotTree(AQ)

figure
plot(AQ)

%% compute lobes
figure; PlotMap3D(AQ, 'lobe');
% Generatlly use the save function when AQ given and returned by function.
% Note that the upper left lobe's lingular is treated seperately. An
% emulation of the middle lobe in the right lung. The Airway segmentation
% must reach all anatomical lobes for this function to be successful.

%% Plot the skeleton inside the segmentation
figure; PlotSegSkel(AQ);

%% traverse the first indexed airway, measure and display results
% this may take >20 minutes
% given the Airway index of interest (use the plottree and plot function to 
% identify branch indices) this function interpolates slices along that
% branch and then attempts to fit an ellipse to the inner lumen, wall peak
% attenuation and outer wall boundaries.
idx = 24;
tic;
CreateAirwayImage(AQ, idx);
FindAirwayBoundariesFWHM(AQ, idx);
disp(toc/60)

% the interpolated slices and fitted ellipses can be viewed with the below
% function. interpolated slice shows 1x1 mm pixels, if the voxel size of
% the original image is greater than this then the interpolated slices will
% be noticibly blurry and results poor.
PlotAirway3(AQ, idx)
OrthoViewAirway(AQ, idx)

% AirQuant doesn't save tapering measurements automatically 
% (only interpolated slices), so we can manually invoke the save method after
% performing measurements on a single branch.
save(AQ)

%% traverse all airways
% this could take a number of hours, results will be saved along the way so
% can be interrupted and rerun at a later time.
tic;
AirwayImageAll(AQ);
% show how long it took
time = toc;
disp(toc/60)

%% compute area of all airways
% This should only take a few minutes and automatically saves once all
% measurements are complete.
FindFWHMall(AQ);
save(AQ)

%% This will tell you which branches have been processed.
% if both the AirwayImageAll and FindFWHMall functions have been run to
% completion successfully then the results of this should show all 1s.
SuccessReport(AQ);

%% See number of airways per generation
figure; AirwayCounts(AQ, 'generation')
figure; AirwayCounts(AQ, 'lobe')

