function AirQuantDir = AirQuantAddPath()
[AirQuantDir,~,~] = fileparts(which('AirQuantAddPath'));
addpath(genpath(AirQuantDir))