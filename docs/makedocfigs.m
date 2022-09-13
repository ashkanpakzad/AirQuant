% Make and save plots for AirQuant's documentation

%% config and load
tubeii = 98;
figdir = 'docs/figs';
if ~exist(figdir,'dir')
    mkdir(figdir)
end

% load
run profiling/CA_base.m

% interp for tube ii
AQnet.tubes(tubeii).MakePatchSlices(AQnet.source, type='source', method='linear');
AQnet.tubes(tubeii).MakePatchSlices(AQnet.seg, type='seg', method='linear');

% fwhm measure for tube ii
num_rays = 180;
ray_interval = 0.2;
AQnet.tubes(tubeii).Measure('AirwayFWHMesl', num_rays, ray_interval);
AQnet.ClassifyLungLobes()

%% make figs
funcs = ["tubeplot", "tubeplot3", "tube_plotspline", "tube_plot3d", ...
    "tube_orthoview", "network_plot", "network_plot3", "network_plot3d",...
    "network_plotspline", "network_orthoview"];

for funs = funcs
    f = figure;
    feval(funs, AQnet, tubeii);
    name = fullfile(figdir,char(funs));
    ext = '.png';
    filename = [name, ext];
    saveas(f,filename)
    close(f)
end


%% tube
% plot
function tubeplot(AQnet, tubeii)
%     AQnet.tubes(tubeii).Plot()
end

% plot3
function tubeplot3(AQnet,tubeii)
    AQnet.tubes(tubeii).Plot3();
end

% plotspline
function tube_plotspline(AQnet,tubeii)
    AQnet.tubes(tubeii).PlotSpline();
end

% plot3D
function tube_plot3d(AQnet,tubeii)
    AQnet.tubes(tubeii).Plot3D();
end

% OrthoView
function tube_orthoview(AQnet,tubeii)
    AQnet.tubes(tubeii).OrthoView();
end


%% network

% plot
function network_plot(AQnet,tubeii)
    AQnet.Plot();
end

% plot3
function network_plot3(AQnet,tubeii)
    AQnet.Plot3();
end

% plot3D
function network_plot3d(AQnet,tubeii)
    AQnet.Plot3D();
end

% plotspline
function network_plotspline(AQnet,tubeii)
    AQnet.PlotSpline(); 
end

function network_orthoview(AQnet,tubeii)
    AQnet.OrthoView();
end

