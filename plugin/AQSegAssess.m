% By Ashkan Pakzad, 2021. ashkanpakzad.github.io
% Parent function to set up and run AQ on several cases based on given
% config file

function AQSegAssess(config)

% set config defaults if they dont exist
checkfield('dataset', 'noset');

checkfield('CTsuf', '_raw.nii.gz');
checkfield('segsuf', '_seg.nii.gz');
checkfield('skelsuf', '_seg_PTKskel.nii.gz');

% Set-up directories
% names input files, give fullpaths if not in matlab path
% add AirQuant library to path
try
    AirQuantDir = AirQuantAddPath();
catch
    warning('Could not find AirQuantAddPath, try running it manually first.')
end

skippedcases = [];

casenames = config.casenames;
disp(['Dataset: ', config.dataset])
for ii = 1:length(casenames)
    casename = char(casenames(ii));
    skip = 0; % do not skip by default
    disp([num2str(ii),' of ', num2str(length(casenames))])
    results_dir = fullfile(AirQuantDir,'results',config.dataset, casename);
    data_dir = fullfile(AirQuantDir,'data', config.dataset);
    logname = [casename, '_log.txt'];
    if ~exist(results_dir, 'dir')
        mkdir(results_dir)
    end
    diary(fullfile(results_dir, logname))
    disp(['Running AirQuant on ', casename]);
    disp(datetime)
    
    % savename is given to automatically save/load results.
    savename = fullfile(results_dir, [casename, '_AQ.mat']);
    
    % check if results already exists
    if exist(savename, 'file')
        warning('AQ object already exists.')
    end
    
    % Get filenames
    CT_name = fullfile(data_dir, [casename, config.CTsuf]);
    seg_name = fullfile(data_dir, [casename, config.segsuf]);
    skel_name = fullfile(data_dir, [casename, config.skelsuf]);
    
    % Check if each datafile exists, if not it throws a warning and skips
    % the current case.
    checkfile(CT_name, 'CT')
    checkfile(seg_name, 'Seg')
    checkfile(skel_name, 'Skel')
    
    if skip == 0
        
        % Load CT data as double
        meta = niftiinfo(CT_name);
        CT = double(niftiread(meta));
        
        % Load Airway segmentation and its skeleton as logicals
        S = logical(niftiread(seg_name));
        skel = logical(niftiread(skel_name));
        
        % Initialise AirQuant
        % Parses CT, segmentation and skeleton to compute airway tree graph
        
        AQ = AirQuant(CT, meta, S, skel, savename);
        
        disp('SEGMENTATION ASSESSOR ONLY');
        
        % Generate initial analysis figures and save
        skelf = figure;
        PlotSegSkel(AQ);
        saveas(skelf, fullfile(results_dir, [casename, '_skel3d.png']))
        saveas(skelf, fullfile(results_dir, [casename, '_skel3d.fig']))
        
        lobef = figure;
        PlotMap3D(AQ, 'lobe');
        saveas(lobef, fullfile(results_dir, [casename, '_lobe3d.png']))
        saveas(lobef, fullfile(results_dir, [casename, '_lobe3d.fig']))
        
        genf = figure;
        PlotMap3D(AQ, 'generation');
        saveas(genf, fullfile(results_dir, [casename, '_gen3d.png']))
        saveas(genf, fullfile(results_dir, [casename, '_gen3d.fig']))
        
        gf = figure;
        plot(AQ);
        saveas(gf, fullfile(results_dir, [casename, '_graph.png']))
        saveas(gf, fullfile(results_dir, [casename, '_graph.fig']))
        
        ACf = figure;
        AirwayCounts(AQ, 'generation', 1)
        saveas(ACf, fullfile(results_dir, [casename, '_AwyCount.png']))
        
        AClf = figure;
        seg_lobegen = AirwayCounts(AQ, 'lobe', 1);
        saveas(AClf, fullfile(results_dir, [casename, '_AwyCountlobe.png']))
        writetable(seg_lobegen, fullfile(results_dir, [casename, '_seg_lobegen.csv']))
        
        disp('SEGMENTATION ASSESSOR ONLY');
      
        % reset
        disp(['Case: ', casename, ' complete.'])
        disp(datetime)
    end
    close all;
    diary off;
end

disp('The following cases were skipped:')
disp(skippedcases)

    function checkfield(field, default)
        % checks if the field for the config struct exists and if not,
        % then sets to the default value provided.
        if ~isfield(config, field)
            config.(field) = default;
        end
    end

    function checkfile(filepath, filetype)
        % checks if the file exists, if not it will throw a warning and
        % inform the script to skip the current case.
        if ~exist(filepath, 'file')
            warning([filetype,' file at path: ', filepath, ' does not exist.'])
            fskip;
        end
    end

    function fskip
        % if informed to skip the current case, it throws a warning if not
        % already for the current case and saves the name to provide a full
        % list at the end.
        if skip == 0
            skippedcases = [skippedcases, string(casename)];
            warning('skipping case...')
        end
        skip = 1;
    end
end

