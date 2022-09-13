% tutorial for running on airways in a clinical chest CT
% By Ashkan Pakzad, 2022. ashkanpakzad.github.io

%% config whole case
AirQuantDir = AirQuantAddPath();
dataset = 'airquant';
casename = 'chestct';
results_dir = fullfile(AirQuantDir,'results', dataset, casename);
savepath = fullfile(AirQuantDir,'profile','chestct_airways_AQ.mat');

CT_name = [casename, '_source.nii.gz'];
seg_name = [casename, '_airway.nii.gz'];
skel_name = [casename, '_airway_PTKskel.nii.gz'];

% CNN methods
sample_sz = 0.5;
modulepath = '/home/ashkan/PhD/awyGAN/';
model_path = '/home/ashkan/PhD/awyGAN/simgan_runs/elated-river-71/checkpoint_9980.tar';

%% Create initial AirQuant Tube Network for Clincal Airways
tic
% Load CT data as double
meta = niftiinfo(CT_name);
source = double(niftiread(meta));

% Load Airway segmentation and its skeleton as logicals
seg = logical(niftiread(seg_name));
skel = logical(niftiread(skel_name));

AQnet = ClinicalAirways(source, meta, seg, skel, fillholes=1, largestCC=1);
AQnet.ClassifyLungLobes()
AQnet.plane_sample_sz = sample_sz;

toc

%% single airways demo
tubeii = 99;

% interp for tube ii
AQnet.tubes(tubeii).MakePatchSlices(AQnet.source, type='source', method='linear', sample_sz=sample_sz);
AQnet.tubes(tubeii).MakePatchSlices(AQnet.seg, type='seg', method='linear', sample_sz=sample_sz);
% measure using model
% AQnet.tubes(tubeii).Measure('AirwayawyGAN', modulepath, model_path);
% figure; AQnet.tubes(tubeiii).plot(smoothing=0.1)
% AQnet.tubes(tubeii).OrthoView()

%% measure using fwhm
num_rays = 60;
ray_interval = 0.2;
AQnet.tubes(tubeiii).Measure('AirwayFWHMesl', num_rays, ray_interval);
figure; AQnet.tubes(tubeiii).OrthoView()
figure; AQnet.tubes(tubeiii).plot()

%% process all airways to make patches
tic
AQnet.MakeTubePatches(method='linear')
toc

%% process CNN measurement
tic
modulepath = '/home/ashkan/PhD/awyGAN/';
model_path = '/home/ashkan/PhD/awyGAN/simgan_runs/elated-river-71/checkpoint_2340.tar';
AQnet.Measure('AirwayawyGAN', modulepath, model_path);
toc
