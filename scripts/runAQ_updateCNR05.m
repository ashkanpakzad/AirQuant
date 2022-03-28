% example config script for batch use of AirQuant by Ashkan Pakzad 2021 
% ashkanpakzad.github.io

%%% add all AirQuant to MATLAB path.
% you may need to provide exact path to this function if it fails.
AirQuantDir = AirQuantAddPath();

% get casenames from config folder
config = [];

% uncomment and specify directory to Input/Output data if not AirQuant/data
% and AirQuant/results
% config.AirQuantDir = 'path/to/dir';

% must be string array, using double dash quotes if specifying cases.
% config.casenames = [""];
config.dataset = 'leuven-ipf';

% suffix for accompanying images
config.CTsuf = '_raw.nii.gz';
config.segsuf= '_seg.nii.gz';
config.skelsuf = '_seg_PTKskel.nii.gz';

% get config file with all casenames in config
config = UpdateResults(config);

%%% pass to AirQuant runner
update_CNR05(config);

