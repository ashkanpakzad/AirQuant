% By Ashkan Pakzad, 2021. ashkanpakzad.github.io
% Parent function to set up and run AQ on several cases based on given
% config file

function skip = wf_clinicalairways_fwhmesl(casename, sourcef, segf, skelf, root_results_dir, gpu)

    % recommended
    num_rays = 180;
    ray_interval = 0.2;

%% run loop
    % init
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
            largestCC=1, spline_sample_sz=0.5, plane_sample_sz=0.5);

        % Generate initial analysis figures and save

        % segskel
        f = figure('WindowStyle', 'Docked'); AQnet.Plot3D(alpha=0.3); hold on; AQnet.Plot3D(type='skel',alpha=1);
        savefig(f,fullfile(results_dir, strcat(casename,"_skel")));
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_skel.png")));
        close(f);

        % plot3
        f = figure('WindowStyle', 'Docked'); AQnet.Plot3D(alpha=0.3); hold on; AQnet.Plot3();
        savefig(f,fullfile(results_dir, strcat(casename,"_plot3")));
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_plot3.png")));
        close(f);

        % splines
        f = figure('WindowStyle', 'Docked'); AQnet.Plot3D(alpha=0.3); hold on; AQnet.PlotSpline();
        savefig(f,fullfile(results_dir, strcat(casename,"_spline")));
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_splines.png")));
        close(f);

        % generation per lobe histogram
        f = figure('WindowStyle', 'Docked'); AQnet.Histogram(label='lobe_gen',region='lobe')
        savefig(f,fullfile(results_dir, strcat(casename,"_hist_lobegen")));
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_hist_lobegen.png")));
        close(f);

        % plot 2d
        f = figure('WindowStyle', 'Docked'); [h, g] = AQnet.Plot(colour='lobe',weightfactor=3);
        try
             AQnet.GraphLobarLayout(h, g)
        catch
        end
        savefig(f,fullfile(results_dir, strcat(casename,"_plot")))
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_plot.png")));
        close(f);

        % lobe 3d + plot 3
        f = figure('WindowStyle', 'Docked'); AQnet.Plot3D(colour='lobe'); hold on; AQnet.Plot3(colour='lobe');
        savefig(f,fullfile(results_dir, strcat(casename,"_plot3d_lobe")));
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_plot3d_lobe.png")));
        close(f);

        % gen 3d
        f = figure('WindowStyle', 'Docked'); AQnet.Plot3D(colour='generation');
        savefig(f,fullfile(results_dir, strcat(casename,"_plot3d_generation")));
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_plot3d_generation.png")));
        close(f);

        % tortuosity
        f = figure('WindowStyle', 'Docked'); AQnet.Histogram(label='tortuosity',region='lobe')
        savefig(f,fullfile(results_dir, strcat(casename,"_hist_tortuosity")));
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_hist_tortuosity.png")));
        close(f);

        % make measurements
        AQnet.MakeTubePatches(method='linear',gpu=gpu)
    
        % make FWHM measurement
        AQnet.Measure('AirwayFWHMesl', num_rays, ray_interval);
    
        % save measures
        save(savename, "AQnet")
        
        % save measurements
        AQnet.ExportCSV(fullfile(results_dir,strcat(casename,"_FWHMesl")))

        % save orthopatches
        AQnet.RunAllTubes('ExportOrthoPatches',fullfile(results_dir,[casename, '_patches']), casename)

        % generate post analysis figures
        % plot 2d - avg D
        f = figure('WindowStyle', 'Docked'); [h, g] = AQnet.Plot(colour='lobe', weight='diameter_mean', ...,
            weightfactor=10, label='lobe_gen');
        try
             AQnet.GraphLobarLayout(h, g)
        catch
        end
        savefig(f,fullfile(results_dir, strcat(casename,"_plot_diameter")))
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_plot_diameter.png")));
        close(f);

        % intertapering
        f = figure('WindowStyle', 'Docked'); AQnet.Histogram(label='intertaper',region='lobe')
        savefig(f,fullfile(results_dir, strcat(casename,"_hist_inter")));
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_hist_inter.png")));
        close(f);

        % diameter
        f = figure('WindowStyle', 'Docked'); AQnet.Histogram(label='diameter_mean',region='lobe')
        savefig(f,fullfile(results_dir, strcat(casename,"_hist_diameter")));
        exportgraphics(f, fullfile(results_dir,strcat(casename,"_hist_diameter.png")));
        close(f);
    
        % reset
        disp(['Case: ', casename, ' complete.'])
        disp(['Total time: ', num2str(toc/60), ' minutes.'])
        disp(datetime)
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