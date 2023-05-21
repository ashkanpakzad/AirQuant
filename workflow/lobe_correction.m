
%% SET PARAMS & LOAD
root_results_dir = 'results';
casename = ['CT079_2'];

results_dir = fullfile(root_results_dir,casename);
savename = fullfile(results_dir, [casename, '_AQnet.mat']);
load(savename)

%% PLOT LOBES
missing_tubes = find(ismissing(string(AQnet.GetTubeValues('lobe', 1))));
disp(missing_tubes)
for ii = missing_tubes
    AQnet.tubes(ii).SetRegion('lobe','T')
end

figure; AQnet.Plot3(colour='lobe');

%% MANIPULATE LOBE SETTINGS
% 'T' , 'B', 'RUL', 'RML', 'RLL', 'LUL', 'LML', 'LLL'
% 
% AQnet.tubes(6).SetRegionDescendants('lobe','LUL')
% AQnet.tubes(6).SetRegion('lobe','B')
% AQnet.tubes(5).SetRegionDescendants('lobe','LUL')
% AQnet.tubes(5).SetRegion('lobe','B')
% AQnet.tubes(11).SetRegionDescendants('lobe','LML')
AQnet.tubes(21).SetRegionDescendants('lobe','LUL')
% AQnet.tubes(12).SetRegionDescendants('lobe','LUL')
% AQnet.tubes(13).SetRegionDescendants('lobe','LUL')
% AQnet.tubes(21).SetRegionDescendants('lobe','LUL')

AQnet.RunAllTubes('SetRegionGeneration','lobe')

%% CHECK LOBES AGAIN
AQnet.RunAllTubes('SetRegionGeneration','lobe')

figure; AQnet.Plot3D(colour='lobe'); hold on; AQnet.Plot3(colour='lobe');
figure; AQnet.Plot(colour='lobe',label='lobe_gen');

%% RESAVE OUTPUTS

AQnet.RunAllTubes('SetRegionGeneration','lobe')
% generation per lobe histogram
f = figure('WindowStyle', 'Docked'); AQnet.Histogram(label='lobe_gen',region='lobe')
savefig(f,fullfile(results_dir, strcat(casename,"_hist_lobegen")));
exportgraphics(f, fullfile(results_dir,strcat(casename,"_hist_lobegen.png")));
close(f);

% plot 2d
f = figure('WindowStyle', 'Docked'); [h, g] = AQnet.Plot(colour='lobe',weightfactor=3);
try
     AQnet.GraphLobarLayout(h, g)
catch
end
savefig(f,fullfile(results_dir, strcat(casename,"_plot")))
exportgraphics(f, fullfile(results_dir,strcat(casename,"_plot.png")));
close(f);

% lobe 3d + plot 3
f = figure('WindowStyle', 'Docked'); AQnet.Plot3D(colour='lobe'); hold on; AQnet.Plot3(colour='lobe');
savefig(f,fullfile(results_dir, strcat(casename,"_plot3d_lobe")));
exportgraphics(f, fullfile(results_dir,strcat(casename,"_plot3d_lobe.png")));
close(f);

% tortuosity
f = figure('WindowStyle', 'Docked'); AQnet.Histogram(label='tortuosity',region='lobe')
savefig(f,fullfile(results_dir, strcat(casename,"_hist_tortuosity")));
exportgraphics(f, fullfile(results_dir,strcat(casename,"_hist_tortuosity.png")));
close(f);

% save OBJECT
save(savename, "AQnet")

% save measurements
AQnet.ExportCSV(fullfile(results_dir,strcat(casename,"_FWHMesl")))

% generate post analysis figures
% plot 2d - avg D
f = figure('WindowStyle', 'Docked'); [h, g] = AQnet.Plot(colour='lobe', weight='diameter_mean', ...,
    weightfactor=10, label='lobe_gen');
try
     AQnet.GraphLobarLayout(h, g)
catch
end
savefig(f,fullfile(results_dir, strcat(casename,"_plot_diameter")))
exportgraphics(f, fullfile(results_dir,strcat(casename,"_plot_diameter.png")));
close(f);

% intertapering
f = figure('WindowStyle', 'Docked'); AQnet.Histogram(label='intertaper',region='lobe')
savefig(f,fullfile(results_dir, strcat(casename,"_hist_inter")));
exportgraphics(f, fullfile(results_dir,strcat(casename,"_hist_inter.png")));
close(f);

% diameter
f = figure('WindowStyle', 'Docked'); AQnet.Histogram(label='diameter_mean',region='lobe')
savefig(f,fullfile(results_dir, strcat(casename,"_hist_diameter")));
exportgraphics(f, fullfile(results_dir,strcat(casename,"_hist_diameter.png")));
close(f);

