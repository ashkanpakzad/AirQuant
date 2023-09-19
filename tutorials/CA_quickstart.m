%% AirQuant Quickstart (Clinical Airways)
% A quickstart tutorial to get going with AirQuant's functionality with
% example data!
tic
% You may need to run `AirQuantAddPath` before running this script to 
% configure MATLAB to run AirQuant.
AirQuantDir = AirQuantAddPath();

%% Download Data
% First download the data!
dataname='chestct';
AQdownload_data(dataname);

%% Load data 
% get data from nifti images into arrays and header variables.

% get CT, segmentation and skeleton images of the case
CT_name = [dataname, '_source.nii.gz'];
seg_name = [dataname, '_airway.nii.gz'];
skel_name = [dataname, '_airway_PTKskel.nii.gz'];

% Load data into relevant formations from nifti images.
% Load the source CT metadata information
meta = niftiinfo(CT_name);
% Load CT data as double
source = niftiread(meta);

% Load Airway segmentation and its skeleton as logicals
seg = logical(niftiread(seg_name));
skel = logical(niftiread(skel_name));


%% Make AirQuant object of data
% There are lots of options here, be sure to check the documentation.
% This will process the skeleton of our image into its components for the AirQuant framework.
% Is the major step where the airway tree is broken down into its constituent branches. 
% Splines are fitted to each branch, making the backbone of our analysis.
AQnet = ClinicalAirways(skel, source=source, header=meta, seg=seg, fillholes=1, largestCC=1, plane_sample_sz=0.5, spline_sample_sz=0.5);
% note: if you are running this in a live matlab script then the progress bar
% will not function properly.

% ClinicalAirways class has its own method for implementing automated lobe classification. We can call this to run on our image.
% The colour option can also be called for any visualisation method. Here
% it is called to show the lobes.
AQnet.ClassifyLungLobes();
%% Basic visualisation
% There are lots of visualisation options. Here we cover some basic options.

figure;
AQnet.Plot3D();

figure; AQnet.Plot3();
figure; AQnet.Plot();
% note that the edge labels correspond to the segment indicies

%% Advanced visualisation
% We can use the lobe region classifications to further enrich our visualisation.
% Checkout the documentation for more advanced visualision use cases.
figure; AQnet.Plot(colour='lobe', weight='generation', weightfactor=10);
title('Edge weighted by generation')
figure; AQnet.Plot3D(colour='lobe');
% `Plot3D` can be resourse demanding on your specs. You maywant to
% skip it.
figure; AQnet.PlotSpline(colour='lobe');

%% Export values to csv
% Finally we can export our values to a csv file for easy analysis in our
% next favourite package!
AQnet.ExportCSV('example.csv');
% print first 10 rows of example.csv
T = readtable('example.csv');
T(1:10,:)
toc
