% By Ashkan Pakzad, 2020. ashkanpakzad.github.io
% This script allows one to set up and run AQ on several cases for one
% script.

%% Set-up directories
% names input files, give fullpaths if not in matlab path
% add AirQuant library to path
AirQuantDir = AirQuantAddPath();


casenames = ["github_demo"];

for ii = 1:length(casenames)
    results_dir = fullfile(AirQuantDir,'results', casenames(ii));
    logname = [char(casenames(ii)), '_log.txt'];
    diary(fullfile(results_dir, logname))
    disp(['Running AirQuant on ', casenames(ii)]);
    if ~exist(results_dir, 'dir')
        mkdir(results_dir)
    end
    
    %% Get filenames
    CT_name = [char(casenames(ii)), '_raw.nii.gz'];
    seg_name = [char(casenames(ii)), '_seg.nii.gz'];
    skel_name = [char(casenames(ii)), '_seg_PTKskel.nii.gz'];

    % Load CT data as double
    meta = niftiinfo(CT_name);
    CT = double(niftiread(meta));

    % Load Airway segmentation and its skeleton as logicals
    S = logical(niftiread(seg_name));
    skel = logical(niftiread(skel_name));

    %% Initialise AirQuant
    % Parses CT, segmentation and skeleton to compute airway tree graph
    % savename is given to automatically save/load results.
    savename = fullfile(results_dir, [char(casenames(ii)), '_AQ.mat']);
    AQ = AirQuant(CT, meta, S, skel, savename);

    %% Generate initial analysis figures and save
    skelf = figure;
    PlotSegSkel(AQ);
    saveas(skelf, fullfile(results_dir, [char(casenames(ii)), '_skel3d.png']))
    saveas(skelf, fullfile(results_dir, [char(casenames(ii)), '_skel3d.fig']))
    
    lobef = figure;
    PlotMap3D(AQ, 'lobe');
    saveas(lobef, fullfile(results_dir, [char(casenames(ii)), '_lobe3d.png']))
    saveas(lobef, fullfile(results_dir, [char(casenames(ii)), '_lobe3d.fig']))
    
    genf = figure;
    PlotMap3D(AQ, 'generation');
    saveas(genf, fullfile(results_dir, [char(casenames(ii)), '_gen3d.png']))
    saveas(lobef, fullfile(results_dir, [char(casenames(ii)), '_lobe3d.fig']))

    gf = figure;
    plot(AQ);
    saveas(genf, fullfile(results_dir, [char(casenames(ii)), '_graph.png']))
    saveas(lobef, fullfile(results_dir, [char(casenames(ii)), '_graph.fig']))
    
    ACf = figure;
    AirwayCounts(AQ, 'generation', 1)
    saveas(ACf, fullfile(results_dir, [char(casenames(ii)), '_AwyCount.png']))
    
    AClf = figure;
    AirwayCounts(AQ, 'lobe', 1)
    saveas(AClf, fullfile(results_dir, [char(casenames(ii)), '_AwyCountlobe.png']))
    
    
    
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
    
    %% Check success
    SuccessReport(AQ);

    close all;
    diary off;
end
