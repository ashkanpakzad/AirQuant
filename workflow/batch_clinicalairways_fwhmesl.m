% By Ashkan Pakzad, 2021. ashkanpakzad.github.io
% Parent function to set up and run AQ on several cases based on given
% config file

function batch_clinicalairways_fwhmesl(varargin)

%% parse input
% dataset dirs
if nargin == 1 || nargin == 2
    dataset = varargin{1};
    source_dir = fullfile(dataset,'source');
    seg_dir = fullfile(dataset,'airway');
    skel_dir = fullfile(dataset,'skel');
elseif nargin == 3 || nargin == 4
    source_dir = varargin{1};
    seg_dir = varargin{2};
    skel_dir = varargin{3};
else
    error('Expected either 1,2,3 or 4 input arguments.')
end

% results dir
if nargin == 2
    results_dir = varargin{2};
elseif nargin == 4
    results_dir = varargin{4};
else % if nargs = 1 or 3 just make in cwd
    results_dir = 'results';
end

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
    
    skip = wf_clinicalairways_fwhmesl(casename, sourcef, segf, skelf, results_dir);
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