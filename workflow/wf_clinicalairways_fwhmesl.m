% Parent function to set up and run AQ on several cases based on given
% config file

function [skip, runinfo] = wf_clinicalairways_fwhmesl(casename, sourcef, segf, skelf, root_results_dir)

    % prune ends of airways
    prune_ends = [2, 2];

    % recommended fwhmesl parameters
    num_rays = 180; 
    ray_interval = 0.2;

%% run loop
    % init
    runinfo = struct('casename',casename);
    skip = 0; % do not skip by default
    disp(['[',casename,'] ','Starting'])
    results_dir = fullfile(root_results_dir,casename);
    mkdir_existok(results_dir)

    % set up log
    logname = fullfile(results_dir,[casename, '_log.txt']);
    diary(logname)
    disp(datetime)
    
    % savename is given to automatically save/load results.
    savename = fullfile(results_dir, [casename, '_AQnet.mat']);
    
    % check if results already exists
    if isfile(savename)
        warning('AQ object already exists.')
    end
    
    % Check if each datafile exists, if not it throws a warning and skips
    % the current case.
    checkfile(sourcef, 'source')
    checkfile(segf, 'segmentation')
    checkfile(skelf, 'skeleton')
    
    if skip == 0
        
        tic; % start timer
        
        % Load CT data as double
        meta = niftiinfo(sourcef);
        CT = double(niftiread(meta));
        
        % Load Airway segmentation and its skeleton as logicals
        S = logical(niftiread(segf));
        skel = logical(niftiread(skelf));
        
        % Initialise AirQuant
        % Parses CT, segmentation and skeleton to compute airway tree graph
        disp(['[',casename,'] ','Init AirQuant.'])
        AQnet = ClinicalAirways(skel, source=CT, header=meta, seg=S,fillholes=1, ...
            largestCC=1, spline_sample_sz=0.5, plane_sample_sz=0.5)

        % save graph
        AQnet.ExportGraph(fullfile(results_dir, [casename, '_graph.csv']));

        %%% Generate initial analysis figures and save
        % segskel
        f = figure(); AQnet.Plot3D(alpha=0.3); hold on; AQnet.Plot3D(type='skel',alpha=1);
        fig_save(f,fullfile(results_dir, strcat(casename,"_skel")));
        close(f);

        % plot3
        f = figure(); AQnet.Plot3D(alpha=0.3); hold on; AQnet.Plot3();
        fig_save(f,fullfile(results_dir, strcat(casename,"_plot3")));
        close(f);

        % splines
        f = figure(); AQnet.Plot3D(alpha=0.3); hold on; AQnet.PlotSpline();
        fig_save(f,fullfile(results_dir, strcat(casename,"_spline")));
        close(f);

        % generation per lobe histogram
        f = figure(); AQnet.Histogram(label='lobe_gen',region='lobe')
        fig_save(f,fullfile(results_dir, strcat(casename,"_hist_lobegen")));
        close(f);

        % plot 2d
        f = figure(); [h, g] = AQnet.Plot(colour='lobe',weightfactor=3);
        try
             AQnet.GraphLobarLayout(h, g)
        catch
        end
        fig_save(f,fullfile(results_dir, strcat(casename,"_plot")))
        close(f);

        % lobe 3d + plot 3
        f = figure(); AQnet.Plot3D(colour='lobe'); hold on; AQnet.Plot3(colour='lobe');
        fig_save(f,fullfile(results_dir, strcat(casename,"_plot3d_lobe")));
        close(f);

        % gen 3d
        f = figure(); AQnet.Plot3D(colour='generation');
        fig_save(f,fullfile(results_dir, strcat(casename,"_plot3d_generation")));
        close(f);

        % tortuosity
        f = figure(); AQnet.Histogram(label='tortuosity',region='lobe')
        fig_save(f,fullfile(results_dir, strcat(casename,"_hist_tortuosity")));
        close(f);

        % set prune
        AQnet.RunAllTubes('SetPruneLength',prune_ends);

        % get patches
        AQnet.MakeTubePatches(method='linear');
    
        % make FWHM measurement
        AQnet.Measure('AirwayFWHMesl', num_rays, ray_interval, 0);
    
        % save measures
        save(savename, "AQnet");
        
        % save measurements
        AQnet.ExportCSV(fullfile(results_dir,strcat(casename,"_FWHMesl")));

        % save orthopatches
        tarpath = fullfile(results_dir,[casename, '_patches.tar']);
        AQnet.RunAllTubes('ExportOrthoPatches',tarpath, casename);
        
        % save grid preview of patches
        f = figure();
        AQnet.grid_preview_measures(4, 4);
        exportgraphics(f,fullfile(results_dir,strcat(casename, '_patch_preview_fwhm.png')));
        close(f);

        % generate post analysis figures
        % plot 2d - avg D
        f = figure(); [h, g] = AQnet.Plot(colour='lobe', weight='diameter_mean', ...,
            weightfactor=10, label='lobe_gen');
        try
             AQnet.GraphLobarLayout(h, g);
        catch
        end
        fig_save(f,fullfile(results_dir, strcat(casename,"_plot_diameter")));
        close(f);

        % intertapering
        f = figure(); AQnet.Histogram(label='intertaper',region='lobe');
        fig_save(f,fullfile(results_dir, strcat(casename,"_hist_inter")));
        close(f);

        % diameter
        f = figure(); AQnet.Histogram(label='diameter_mean',region='lobe');
        fig_save(f,fullfile(results_dir, strcat(casename,"_hist_diameter")));
        close(f);
    
        % reset
        disp(['Case: ', casename, ' complete.'])
        disp(['Total time: ', num2str(toc/60), ' minutes.'])
        disp(datetime)

        % save info to summary table
        % total number of tubes
        runinfo.n = length(AQnet.tubes);
        % lengths
        lengths = cell2mat(AirQuant.list_property({AQnet.tubes.stats},'arclength'));
        runinfo.arclen_mm = sum(lengths);        
        % lumen volume
        vols = cell2mat(AirQuant.list_property({AQnet.tubes.stats},'volume'))/1e6;
        vols = vols(1,:);
        runinfo.lumen_vol_l = sum(vols,'omitnan');
        % max lobe generation
        gens = cell2mat(AirQuant.list_property({AQnet.tubes.region},'lobe_gen'));
        runinfo.maxlobegen = max(gens);
        runinfo.runtime_m = toc/60;
        % breakdown by lobe
        try
            lobes = AirQuant.list_property({AQnet.tubes.region},'lobe');
            lobecodes = {'RUL','RML','RLL','LUL','LML','LLL'};
            for i = 1:length(lobecodes)
                idx = strcmp(lobes,lobecodes{i});
                runinfo.([lobecodes{i},'_n']) = sum(idx);
                runinfo.([lobecodes{i},'_len']) = sum(lengths(idx));
                runinfo.([lobecodes{i},'_vol']) = sum(vols(idx),'omitnan');
                runinfo.([lobecodes{i}, '_maxlobegen']) = max(gens(idx));
            end
        catch
            warning([casname,' lobe info not available'])
        end
        % append runinfo to csv
        runinfo_T = struct2table(runinfo);
        writetable(runinfo_T, fullfile(root_results_dir, 'summary.csv'), 'WriteMode', 'append');
    end



    close all;
    diary off;

    function checkfile(filepath, filetype)
        % checks if the file exists, if not it will throw a warning and
        % inform the script to skip the current case.
        if ~exist(filepath, 'file')
            warning([filetype,' file at path: ', filepath, ' does not exist.'])
        end
    end

end