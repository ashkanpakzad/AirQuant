%% parse input
top_results_dir = 'results';

% get list of casenames from source dir
alldir = dir(fullfile(top_results_dir));
% remove any that begin with '.' (i.e. '.', '..', and '.DS_Store')
alldir = alldir(~startsWith({alldir.name}, '.'));
casenames = string({(alldir.name)});

for ii = 28:length(casenames)
    casename = casenames(ii);
    disp(casename)
    disp([num2str(ii), ' of ', num2str(length(casenames))])
    results_dir = fullfile(top_results_dir, casename);
    % load data
    load(fullfile(results_dir, casename + '_AQnet.mat'));
    AQnet.ExportGraph(fullfile(results_dir, casename + '_graph.csv'));
    clear AQnet
end