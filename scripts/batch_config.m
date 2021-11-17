% example config script for batch use of AirQuant by Ashkan Pakzad 2021 
% ashkanpakzad.github.io

%%% add all AirQuant to MATLAB path.
% you may need to provide exact path to this function if it fails.
AirQuantDir = AirQuantAddPath(); 

% get casenames from config folder
config = [];
% must be string array, using double dash quotes.
config.casenames = ["example"];
config.dataset = 'example';

% suffix for accompanying images
config.CTsuf = '_raw.nii.gz';
config.segsuf= '_seg.nii.gz';
config.skelsuf = '_seg_PTKskel.nii.gz';

%%% pass to AirQuant runner
runAQ(config);

