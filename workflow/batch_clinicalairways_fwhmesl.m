% By Ashkan Pakzad, 2021. ashkanpakzad.github.io
% Parent function to set up and run AQ on several cases based on given
% config file

function batch_clinicalairways_fwhmesl(source_dir, seg_dir, skel_dir, gpu)

%% parse input
results_dir = 'results';

% get list of casenames from source dir
alldir = dir(fullfile(source_dir));
alldir = alldir(~[alldir.isdir]);
raw_names = string({(alldir.name)});
casenames = strrep(raw_names,'.nii.gz','') ;


%% Set-up directories
% names input files, give fullpaths if not in matlab path
% add AirQuant library to path
try
    AirQuantAddPath();
catch
    warning('Could not find AirQuantAddPath, try running it manually first.')
end

mkdir_existok(results_dir)

%% run loop
skippedcases = [];

for ii = 1:length(casenames)
    % init
    casename = char(casenames(ii));
    skip = 0; % do not skip by default
    disp([num2str(ii),' of ', num2str(length(casenames))])
    
    % Get filenames
    filesuffix = '.nii.gz';
    sourcef = fullfile(source_dir,[casename,filesuffix]);
    segf = fullfile(seg_dir,[casename,filesuffix]);
    skelf = fullfile(skel_dir,[casename,filesuffix]);
    
    skip = wf_clinicalairways_fwhmesl(casename, sourcef, segf, skelf, results_dir, gpu);
    if skip
        fskip;
    end
end

disp('The following cases were skipped:')
disp(skippedcases)

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