% profiling - loading base
% add AirQuant library to path
tic
AirQuantDir = AirQuantAddPath();
dataset = 'example';
casename = 'github_demo';
results_dir = fullfile(AirQuantDir,'results', dataset, casename);
savepath = fullfile(AirQuantDir,'profile','profile_AQ.mat');
% By Ashkan Pakzad, 2022. ashkanpakzad.github.io

%% load and init
CT_name = [casename, '_raw.nii.gz'];
seg_name = [casename, '_seg.nii.gz'];
skel_name = [casename, '_seg_PTKskel.nii.gz'];

% % Load CT data as double
meta = niftiinfo(CT_name);
CT = double(niftiread(meta));

% Load Airway segmentation and its skeleton as logicals
S = logical(niftiread(seg_name));
skel = logical(niftiread(skel_name));

AQnet = TubeNetwork(CT, meta, S, skel);
toc