% profiling - loading base
% add AirQuant library to path
tic
AirQuantDir = AirQuantAddPath();
casename = 'chestct';
dataset = AQdownload_data(casename);
AirQuantAddPath();
savepath = fullfile(AirQuantDir,'profile','profile_AQ.mat');
% By Ashkan Pakzad, 2022. ashkanpakzad.github.io

%% load and init
CT_name = [casename, '_source.nii.gz'];
seg_name = [casename, '_airways.nii.gz'];
skel_name = [casename, '_airways_PTKskel.nii.gz'];

% % Load CT data as double
meta = niftiinfo(CT_name);
CT = double(niftiread(meta));

% Load Airway segmentation and its skeleton as logicals
S = logical(niftiread(seg_name));
skel = logical(niftiread(skel_name));

AQnet = TubeNetwork(skel, seg=seg, source=source, header=meta, fillholes=1, largestCC=1);
toc