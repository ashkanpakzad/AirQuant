% By Ashkan Pakzad, 2020. ashkanpakzad.github.io

%% Set-up directories
% names input files, give fullpaths if not in matlab path
% add AirQuant library to path
AirQuantDir = AirQuantAddPath();
results_dir = fullfile(AirQuantDir,'results');
casename = 'github_demo';

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
AQ = ComputeAirwayLobes(AQ);
figure; PlotMap3D(AQ, 'Lobe');
% Generatlly use the save function when AQ given and returned by function.
save(AQ) 
% MATLAB's colourcoding can be a little buggy when viewing results.
% Note that the upper left lobe's lingular is treated seperately. An
% emulation of the middle lobe in the right lung. The Airway segmentation
% must reach all anatomical lobes for this function to be successful.

%% Plot the skeleton inside the segmentation
figure
patch(isosurface(S),'EdgeColor', 'none','FaceAlpha',0.3);
hold on
isosurface(skel)
axis vis3d

%% traverse the first indexed airway, measure and display results
% this may take >20 minutes
% given the Airway index of interest (use the plottree and plot function to 
% identify branch indices) this function interpolates slices along that
% branch and then attempts to fit an ellipse to the inner lumen, wall peak
% attenuation and outer wall boundaries.
idx = 26;
tic;
AQ = CreateAirwayImage(AQ, idx);
AQ = FindAirwayBoundariesFWHM(AQ, idx);
disp(toc/60)

% the interpolated slices and fitted ellipses can be viewed with the below
% function. interpolated slice shows 1x1 mm pixels, if the voxel size of
% the original image is greater than this then the interpolated slices will
% be noticibly blurry and results poor.
PlotAirway3(AQ, idx)

% AirQuant doesn't save tapering measurements automatically 
% (only interpolated slices), so we can manually invoke the save method after
% performing measurements on a single branch.
save(AQ)

%% traverse all airways
% this could take a number of hours, results will be saved along the way so
% can be interrupted and rerun at a later time.
tic;
AQ = AirwayImageAll(AQ);
% show how long it took
time = toc;
disp(toc/60)

%% compute area of all airways
% This should only take a few minutes and automatically saves once all
% measurements are complete.
AQ = FindFWHMall(AQ);
save(AQ)

%% This will tell you which branches have been processed.
% if both the AirwayImageAll and FindFWHMall functions have been run to
% completion successfully then the results of this should show all 1s.
report = debuggingreport(AQ);

%% Compute taper gradient path of all terminal branches
% Table with all carina-terminal branch data.
% running ComputeAirwayLobes() before this function will group output by lobe.
AllTaperResults = ComputeTaperAll(AQ);

%% Display taper results of a single gradient path
% plot taper grad results from carina to a given end-node.
figure
PlotTaperResults(AQ, AllTaperResults.terminalnode(1), 'inner')
figure
PlotTaperResults(AQ, AllTaperResults.terminalnode(1))

%% Construct single taper gradient path
% this demonstrates more hands on functions to process data as desired.
% list of terminal node, refer to airway graph.
terminalnodelist = ListTerminalNodes(AQ);
% get branch data for carina to terminal node for given node.
% note a positive taper gradient shows that the branch is tapering.
[logtaperrate, cum_arclength, cum_area, path] = ConstructTaperPath(AQ, terminalnodelist(1)); 

%% Display boxplot of tapervalues by lobe
figure
TaperBoxPlot(AQ, 'inner')
figure
TaperBoxPlot(AQ)
