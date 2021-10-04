% By Ashkan Pakzad (ashkanpakzad.github.io) June 2021 
function config = RunWholeDataset(config)
% user can provide config structure with dataset name and this function
% will populate the 'casenames' field with all the cases in that dataset
% based on the source CT images available.

% get list of file names in dataset folder
AirQuantDir = AirQuantAddPath(); 
alldir = dir(fullfile(AirQuantDir,'data',config.dataset));
alldir = alldir(~[alldir.isdir]);

% identify CTs in dataset
rawfiles_log = false(length(alldir),1);
for ii = 1:length(alldir)
   rawfiles_log(ii) = contains(alldir(ii).name, config.CTsuf);
end

% get list of casenames
raw_names = string({(alldir(rawfiles_log).name)});

% save casenames
config.casenames = strrep(raw_names,config.CTsuf,'') ;

end