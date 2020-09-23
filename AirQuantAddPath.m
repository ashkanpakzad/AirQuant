function AirQuantDir = AirQuantAddPath()
[AirQuantDir,~,~] = fileparts(which('AirQuantAddPath'));
addpath(genpath(AirQuantDir))

results_dir = fullfile(AirQuantDir, 'results');

if ~isfolder(results_dir)
    mkdir(results_dir)
end
end
