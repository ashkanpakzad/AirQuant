% Script by Ashkan Pakzad to compare intra/intertapering across lobes

%% Load data

AirQuantDir = AirQuantAddPath();
results_dir = fullfile(AirQuantDir,'results');
casename = 'N1';

load('N1_segmenttaper_0prune.mat', 'SegmentTaperResults')

%% Boxplot

figure
logtaperdata = [N1_0prune.inner_intra];
typetxt = 'Inner lumen intratapering';
datalabels = [N1_0prune.lobe];

boxplot(logtaperdata, datalabels)
xlabel('Lobe')
ylabel('Intratapering (%)')
title(typetxt)

figure
logtaperdata = [N1_0prune.inner_inter];
typetxt = 'Inner lumen intertapering';
datalabels = [N1_0prune.lobe];

boxplot(logtaperdata, datalabels)
xlabel('Lobe')
ylabel('Intertapering (%)')
title(typetxt)

