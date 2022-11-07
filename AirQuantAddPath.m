function AirQuantDir = AirQuantAddPath(AirQuantDirIn)
% Add AirQuant's directories to the MATLAB path so that its libraries can be used.
% Optional input argument to the directory that user may prefer to store
% data and save outputs.
% AirQuantDir out argument is the directory configured for Input/Output
% See <a href = "https://airquant.readthedocs.io/">AirQuant Documentation</a> for help.'

disp('See <a href = "https://airquant.readthedocs.io/">AirQuant Documentation</a> for help.')

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
