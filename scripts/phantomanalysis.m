% By Ashkan Pakzad October 2020 ashkanpakzad.github.io

% Phantom analysis with AirQuant.

%% Set-up directories
% names input files, give fullpaths if not in matlab path
% add AirQuant library to path
AirQuantDir = AirQuantAddPath();
results_dir = fullfile(AirQuantDir,'results');
recontype = 'Lung'; % phantom reconstruction kernel. Lung or Body.
phantomtype = 'constant'; % taper or constant

samplesizes = 0.1:0.1:1;

% fixed points for taper tubes
% I = [   335, 214, 179; 335, 215, 124;
%         294, 201, 179; 295, 201, 124;
%         350, 174, 179; 349, 176, 124;
%         375, 227, 179; 372, 226, 124;
%         323, 251, 179; 323, 249, 124;   ];

% fixed points for constant tubes
I = [   243, 362, 179; 246, 362, 124;
    213, 360, 179; 215, 360, 124;
    184, 359, 179; 187, 358, 124;
    215, 333, 179; 217, 333, 124;
    211, 387, 179; 213, 387, 124;   ];

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
skel = logical(SkelPTK_windex(S, fixedpoints));

%% Initialise AirQuant
% Parses CT, segmentation and skeleton to compute airway tree graph
% savename is given to automatically save/load results.
for ii = 1:length(samplesizes)
    plane_sample_sz = samplesizes(ii);
    text = sprintf('%1.1f',plane_sample_sz);
    extra = ['planesample_',text(1),'_',text(3)];
    
    casename = ['Phantom_',recontype,'_',phantomtype,'_',extra];
    
    savename = fullfile(results_dir, [casename, '_AQ.mat']);
    AQ = AirQuantPhantom(CT, meta, S, skel, savename);
    AQ.plane_sample_sz = plane_sample_sz;
    
    %% Display Airway Graph in 3D and 2D
    % figure
    % PlotTree(AQ)
    %
    % figure
    % plot(AQ,'gen')
    %
    % Plot the skeleton inside the segmentation
    % skel = logical(Skeletonise_via_PTK_with_index(S, fixedpoints));
    % figure
    % patch(isosurface(S),'EdgeColor', 'none','FaceAlpha',0.3);
    % hold on
    % isosurface(skel)
    % axis vis3d
    % ax = gca;
    % ax.XDir = 'reverse';
    
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
end
