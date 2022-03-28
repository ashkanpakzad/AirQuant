% By Ashkan Pakzad, 2022. ashkanpakzad.github.io
% rerun existing case and evaluate using CNR, save new taperresults

function update_prune25(config)

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
    
    if skip == 0
        
        AQ = AirQuant(savename);

        obj.MeasureMode = 'FWHM';
        
        % save segment taper analysis to csv
        SegmentTaperResults = SegmentTaperAll(AQ, [2.5 2.5]);
        writetable(SegmentTaperResults, fullfile(results_dir, [casename, '_SegmentTaper.csv']));

        % reset
        disp(['Updated Case: ', casename, ' to CNReval'])
        disp(datetime)
        
        clear vars AQ
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


