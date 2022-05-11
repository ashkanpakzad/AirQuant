% profiling - loading base 
% add AirQuant library to path
% CLINICAL AIRWAYS
tic
AirQuantDir = AirQuantAddPath();
dataset = 'airquant';
casename = 'chestct';
results_dir = fullfile(AirQuantDir,'results', dataset, casename);
savepath = fullfile(AirQuantDir,'profile','chestct_airways_AQ.mat');
% By Ashkan Pakzad, 2022. ashkanpakzad.github.io

%% load and init
CT_name = [casename, '_source.nii.gz'];
seg_name = [casename, '_airway.nii.gz'];
skel_name = [casename, '_airway_PTKskel.nii.gz'];

% % Load CT data as double
meta = niftiinfo(CT_name);
CT = double(niftiread(meta));

% Load Airway segmentation and its skeleton as logicals
S = logical(niftiread(seg_name));
skel = logical(niftiread(skel_name));

AQnet = ClinicalAirways(CT, meta, S, skel, fillholes=1, largestCC=1);
toc