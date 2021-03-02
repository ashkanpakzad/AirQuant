% example config script for AirQuant by Ashkan Pakzad 2021 
% ashkanpakzad.github.io

%%% add all AirQuant to MATLAB path.
% you may need to provide exact path to this function if it fails.
AirQuantDir = AirQuantAddPath(); 

%%% set up the config structure
config = [];
% must be string array, using double dash quotes.
config.casenames = ["github_demo"]; % required
% results will be saved in AirQuant/results/example/
config.dataset = 'example';

% suffix of each file (these are defaults if not provided)
config.CTsuf = '_raw.nii.gz';
config.segsuf= '_seg.nii.gz';
config.skelsuf = '_seg_PTKskel.nii.gz';

%%% pass to AirQuant runner
runAQ(config);

