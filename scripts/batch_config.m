% example config script for batch use of AirQuant by Ashkan Pakzad 2021 
% ashkanpakzad.github.io

%%% add all AirQuant to MATLAB path.
% you may need to provide exact path to this function if it fails.
AirQuantDir = AirQuantAddPath(); 

% get casenames from config folder
config = [];
% must be string array, using double dash quotes.
config.casenames = []; 
config.dataset = 'test';

% suffix for accompanying images
config.CTsuf = '_raw.nii.gz';
config.segsuf= '_seg.nii.gz';
config.skelsuf = '_seg_PTKskel.nii.gz';

% get config file with all casenames in config
config = RunWholeDataset(config);

%%% pass to AirQuant runner
AQSegAssess(config);

