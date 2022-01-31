function casenames = GetDatasetCasenames(config)

suffix = '_raw.nii.gz';

AirQuantDir = AirQuantAddPath();
ds_path = fullfile(AirQuantDir, 'data', config.dataset);

casedir = dir(ds_path);

% remove . and ..
casedir = casedir(~ismember({casedir.name},{'.','..'}));

% reduce to _raw.nii.gz only
casedir_raw = casedir(endsWith({casedir.name},suffix));

% output names, remove suffix
ccasenames = {casedir_raw.name};
ccasenames = strrep(ccasenames,suffix,'') ;
casenames = string(ccasenames);
end