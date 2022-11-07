%% AirQuant FWHM (Clinical Airways)
% This tutorial assumes you have run the quickstart tutorial and have built
% the `ClinicalAirways` object saved as `AQnet`.
%
% Now that the backbone of the AirQuant analysis has been set up, we can work through subsequent steps of interpolating the source image for airway measurement. 
% we opt to only demonstrate this on one airway, though in reality we want
% this to be done on every airway. We pick one airway id.
tic
% You may need to run `AirQuantAddPath` before running this script to 
% configure MATLAB to run AirQuant.
AirQuantDir = AirQuantAddPath();

%% Single airway/tube patch interpolation and FWHMesl
% we can access any 'tube' (airway) by its index, this allows us to access
% the methods of tubes.
atube = AQnet.tubes(10);

% We run the method make the airway patch slices of this airway on the CT.
usegpu = 0;
atube.MakePatchSlices(AQnet.source, type='source', method='linear', gpu=usegpu);

% set up parameters of FWHMesl
% With our airway patches interpolated we now run the measurement of each patch. We a very straightforward method, FWHMesl.
num_rays = 60;
ray_interval = 0.2;
% segmentation interpolation is also required for the FWHMesl method
atube.MakePatchSlices(AQnet.seg, type='seg', method='linear', gpu=usegpu);
atube.Measure('AirwayFWHMesl', num_rays, ray_interval);

%   We can interactively visualise the CT image along the airway patches with its estimated ellipse fitting.
% (opens externally)
figure; atube.OrthoView();

%% Basic visualisation
% We can plot any measurement along the airway against its arc-length
figure; atube.plot(Y='diameters');
% though this is the default options for X and Y, we will be explicit.

%% Export individual airway/tube
atube.ExportCSV('example.csv');
% print first 10 rows of example.csv
T = readtable('example.csv');
T(1:10,:)


%% process all airways to make patches
% We can process all airways now. This can take somewhere between 1-4 hours to run. 
% By default it will try to run on GPU if available which can speed things up.
% uncomment the lines below to run.

% AQnet.MakeTubePatches(method='linear', gpu=usegpu)
% AQnet.Measure('AirwayFWHMesl', num_rays, ray_interval);

% We can now visualise the characteristics derived from diameter measurements.
% figure; AQnet.Plot(colour='lobe', weight='meandiameter', weightfactor='10')

toc
