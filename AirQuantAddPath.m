function varargout = AirQuantAddPath(AirQuantDirIn)
% Add AirQuant's directories to the MATLAB path so that its libraries can be used.
% Optional input argument to the directory that user may prefer to store
% data and save outputs.
% AirQuantDir out argument is the directory configured for Input/Output
% See <a href = "https://airquant.readthedocs.io/">AirQuant Documentation</a> for help.'

disp('See <a href = "https://airquant.readthedocs.io/">AirQuant Documentation</a> for help.')

% check if required toolboxes are installed
required_packages = {'Signal Processing Toolbox', 'Image Processing Toolbox', ...
    'Statistics and Machine Learning Toolbox', 'Curve Fitting Toolbox',...
    'Parallel Computing Toolbox'};
for apackage = required_packages
    if contains(AQstruct2array(ver), apackage) == 0
        warning(strcat(apackage{1},' is not installed, some features may not work properly.'))
    end
end


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

if nargout > 0
    varargout{1} = AirQuantDir;
end
% configure expected paths for AirQuant
results_dir = fullfile(AirQuantDir, 'results');

% make results directory
if ~isfolder(results_dir)
    mkdir(results_dir)
end

function a = AQstruct2array(s)
% AQSTRUCT2ARRAY Convert structure with doubles to an array. 
% Incase struct2array is not available.
% 
% From
% https://uk.mathworks.com/matlabcentral/answers/1717910-was-struct2array-removed-from-matlab
%

% Convert structure to cell
c = struct2cell(s);
% Construct an array
a = [c{:}];

end
end
