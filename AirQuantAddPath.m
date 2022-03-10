function AirQuantDir = AirQuantAddPath(AirQuantDirIn)
% Add AirQuant's directories to the MATLAB path so that its libraries can be used.
% Optional input argument to the directory that user may prefer to store
% data and save outputs.
% AirQuantDir out argument is the directory configured for Input/Output

% example expected structure:
% AirQuantDir/
    % data/
    %   exampledataset/
    %       examplecase1/
    %           examplecase1_raw.nii.gz
    %           examplecase1_seg.nii.gz
    %           examplecase1_seg_PTKskel.nii.gz
    %           examplecase2_raw.nii.gz
    %           examplecase2_seg.nii.gz
    %           examplecase2_seg_PTKskel.nii.gz
    % results/
    %   exampledataset/
    %       examplecase1/
    %           examplecase1_AQ.mat
    %       examplecase2/
    %          examplecase2_AQ.mat


% get path to airquant library directory and add to matlab path
[AirQuantDirLib,~,~] = fileparts(which('AirQuantAddPath'));
addpath(genpath(AirQuantDirLib))

if nargin < 1
    % use default path if no argument
    AirQuantDir = AirQuantDirLib;
else
    % if data stored outside of airquant directory, add its path to matlab
    addpath(genpath(AirQuantDirIn))
    AirQuantDir = AirQuantDirIn;
end

% configure expected paths for AirQuant
results_dir = fullfile(AirQuantDir, 'results');

% make results directory
if ~isfolder(results_dir)
    mkdir(results_dir)
end
end
