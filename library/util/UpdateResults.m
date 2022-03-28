% By Ashkan Pakzad (ashkanpakzad.github.io) June 2021 
function config = UpdateResults(config)
% user can provide config structure with dataset name and this function
% will populate the 'casenames' field with all the cases in that dataset
% based on the results dir

% get list of file names in dataset folder
AirQuantDir = AirQuantAddPath(); 
files = dir(fullfile(AirQuantDir,'results',config.dataset));

% remove . and .. dir
folderNames = {files([files.isdir]).name};
folderNames = folderNames(~ismember(folderNames ,{'.','..'}));

% save casenames
config.casenames = convertCharsToStrings(folderNames) ;

end