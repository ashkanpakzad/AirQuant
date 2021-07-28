% export analysis tables for each case of a dataset
% By Ashkan Pakzad (ashkanpakzad.github.io) July 2021

% dataset details
AirQuantDir = AirQuantAddPath();

config = [];
config.dataset = 'test';
config.casenames = [];

% get casenames based on dataset files in 'data' folder
casenames = GetDatasetCasenames(config);

% execute function to merge tables
table_data = LoadTaperTables(config); 

% write tables
writetable(table_data, fullfile(AirQuantDir, 'results', config.dataset,'allairways.csv'))

function fulltable =  LoadTaperTables(config)
AirQuantDir = AirQuantAddPath();

dataset = config.dataset;
casenames = config.casenames;

fulltable = [];

for ii = 1:length(casenames)
    % set up paths and log for case
    casename = char(casenames(ii));
    results_dir = fullfile(AirQuantDir, 'results', dataset, casename);

    %% Load taper table object
    newtable = readtable(fullfile(results_dir, [casename, '_SegmentTaper.csv']));
    
    casenamestab = table(repmat(string(casename),height(newtable),1));
    % add casename col
    newtable = [casenamestab newtable];
    % concatenate underneath
    fulltable = [fulltable; newtable];
end

fulltable.Properties.VariableNames{'Var1'} = 'case';
end