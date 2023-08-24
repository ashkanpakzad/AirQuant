% Parent function to set up and run AQ on several cases based on given
% config file

function batch_clinicalvessels_fwhmesl(source_dir, seg_dir, skel_dir, overwrite)

if nargin < 1
    % default behaviour is not to overwrite.
    source_dir = 'source';
end

if nargin < 2
    seg_dir = 'vessels';
end

if nargin < 3
    skel_dir = 'vessels_skel';
end

if nargin < 4
    % default behaviour is not to overwrite.
    overwrite = 0;
end

%% parse input
results_dir = 'results';

% get list of casenames from source dir
alldir = dir(fullfile(source_dir));
alldir = alldir(~[alldir.isdir]);
raw_names = string({(alldir.name)});
casenames = strrep(raw_names,'.nii.gz','');


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

    % check if it already exists
    case_dir = fullfile(results_dir,casename);
    if overwrite == 0 && isfolder(case_dir)
        disp(['[',casename,'] already exists. overwrite == false, therefore skipping.'])
        fskip;
        continue
    end
    
    % Get filenames
    filesuffix = '.nii.gz';
    sourcef = fullfile(source_dir,[casename,filesuffix]);
    segf = fullfile(seg_dir,[casename,filesuffix]);
    skelf = fullfile(skel_dir,[casename,filesuffix]);
    try
    skip = wf_clinicalvessels_fwhmesl(casename, sourcef, segf, skelf, results_dir);
    catch e %e is an MException struct
        fprintf(2,'error identifier :\n%s\n',e.identifier);
        fprintf(2,'There was an error! The message was:\n%s\n',e.message);
        skip = 1;
    end

    if skip == 1
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