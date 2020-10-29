% By Ashkan Pakzad October 2020 ashkanpakzad.github.io

% Phantom analysis with AirQuant.

%% Set-up directories
% names input files, give fullpaths if not in matlab path
% add AirQuant library to path
AirQuantDir = AirQuantAddPath();
results_dir = fullfile(AirQuantDir,'results');
recontype = 'Lung'; % phantom reconstruction kernel. Lung or Body.
phantomtype = 'taper';
casename = ['Phantom_',recontype,'_',phantomtype];

% fixed points for taper tubes
I = [   335, 214, 179; 335, 215, 124; 
        294, 201, 179; 295, 201, 124; 
        350, 174, 179; 349, 176, 124; 
        375, 227, 179; 372, 226, 124;
        323, 251, 179; 323, 249, 124;   ];

%% Get filenames

CT_name = ['Phantom_',recontype, '_raw.nii.gz'];
seg_name = ['Phantom_',phantomtype, '_seg.nii.gz'];

% Load CT data as double
meta = niftiinfo(CT_name);
CT = double(niftiread(meta));

% Load Airway segmentation and its skeleton as logicals
S = logical(niftiread(seg_name));

% convert skel fixed indices to linear idx
fixedpoints = sub2ind(size(CT), I(:,1), I(:,2), I(:,3));
% get skel using basic PTK method
skel = logical(Skeletonise_via_PTK_with_index(S, fixedpoints));

%% Initialise AirQuant
% Parses CT, segmentation and skeleton to compute airway tree graph
% savename is given to automatically save/load results.
savename = fullfile(results_dir, [casename, '_AQ.mat']);
AQ = AirQuantPhantom(CT, meta, S, skel, savename);

%% Display Airway Graph in 3D and 2D
figure
PlotTree(AQ)

%%
figure
plot(AQ,'gen')

%% Plot the skeleton inside the segmentation
figure
patch(isosurface(S),'EdgeColor', 'none','FaceAlpha',0.3);
hold on
isosurface(skel)
axis vis3d
ax = gca;
ax.XDir = 'reverse';

%% traverse the first indexed airway, measure and display results
% this may take >20 minutes
% given the Airway index of interest (use the plottree and plot function to 
% identify branch indices) this function interpolates slices along that
% branch and then attempts to fit an ellipse to the inner lumen, wall peak
% attenuation and outer wall boundaries.
idx = 1;
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
%%
PlotAirway3(AQ, 2)

%%
idx = 2;
tic;
AQ = CreateAirwayImage(AQ, idx);
AQ = FindAirwayBoundariesFWHM(AQ, idx);
disp(toc/60)

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

