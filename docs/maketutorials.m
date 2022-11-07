% Make and save tutorials for AirQuant's documentation

%% set paths
AirQuantDir = AirQuantAddPath();
tutdir = fullfile(AirQuantDir,'tutorials');
outpath = fullfile(AirQuantDir,'docs','_static');
% must be run in specific order
tutfiles = ["CA_quickstart.m","CA_FWHMesl.m"];

options.format = 'html';
options.showCode = true;
options.outputDir = outpath;

%% make tutorials
for ii = 1:length(tutfiles)
    filename = fullfile(tutdir,tutfiles(ii));
    publish(char(filename),options)
    close all
end
