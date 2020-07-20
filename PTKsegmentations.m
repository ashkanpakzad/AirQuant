% By Ashkan Pakzad on 19 July ashkanpakzad.github.io

% Get Airway segmentation, Airway 'processed' skeleton, lungs, lobes via
% PTK. Also converts original dicom to nifti format.

%% Set up
casename = 'M1';
savepath = '/Users/apakz/PhD/normal_data/M1/PTKoutput';
overwritelog = 1; % If log should be reset if already exists (0 to append)

% dicom directory containing only the series to segment
dcmdir = '/Users/apakz/PhD/normal_data/M1/dcm';

% path to PTK source 
PTK_source = '/Users/apakz/PhD/publicPTK';


%% Check for folder and log files
initial_path = pwd;
% create savepath dir if not exists
if ~isfolder(savepath)
    mkdir(savepath)
end

% delete log file if already exists
if overwritelog == 1 && isfile(fullfile(savepath,'corereporting.log'))
    cd(savepath)
    delete corereporting.log
    cd(initial_path)
end

%% Call PTK and load image
% add PTK to path
cd(PTK_source)
PTKAddPaths;
cd(initial_path)

% set up reporting, add to log
reporting = CoreReporting([], true, fullfile(savepath,'corereporting.log'));
ptk_main = PTKMain(reporting);
reporting.Log('Data processed using script by Ashkan Pakzad ashkanpakzad.github.io')
reporting.Log('PulmonaryToolKit(PTK) by Tom Doel github.com/tomdoel/pulmonarytoolkit')
reporting.Log('PTK to process chest CT to return Lung, Lobe, Airway segmentations and Airway centreline')
reporting.Log(['casename: ', casename])
reporting.LogVerbose(['PTK source: ', PTK_source])
reporting.LogVerbose(['savepath: ', savepath])

% load data
reporting.LogVerbose('Loading data')
dataset = ptk_main.Load(dcmdir);

%% Execute PTK modules
% get image volume
reporting.LogVerbose('Starting PTKOriginalImage plugin')
image = dataset.GetResult('PTKOriginalImage');
reporting.LogVerbose('Complete PTKOriginalImage plugin')

% lung segmentation
reporting.LogVerbose('Starting PTKLeftAndRightLungs plugin')
lungs = dataset.GetResult('PTKLeftAndRightLungs');
lungs.Title = [casename, ' / ', 'Segmentation', ' / ', 'PTKLeftAndRightLungs'];
reporting.LogVerbose('Complete PTKLeftAndRightLungs plugin')

% Airway segmentation
reporting.LogVerbose('Starting PTKAirways plugin')
airway_results = dataset.GetResult('PTKAirways');
reporting.LogVerbose('Complete PTKAirways plugin')

% Airway centreline extraction
reporting.LogVerbose('Starting PTKAirwayCentreline plugin')
skeleton_results = dataset.GetResult('PTKAirwayCentreline');
reporting.LogVerbose('Complete PTKAirwayCentreline plugin')

% Lobe segmentation
reporting.LogVerbose('Starting PTKLobes plugin')
lobes = dataset.GetResult('PTKLobes');
lobes.Title = [casename, ' / ', 'Segmentation', ' / ', 'PTKLobes'];
reporting.LogVerbose('Complete PTKLobes plugin')

%% Process airway results into binary matrix as PTKimage
reporting.LogVerbose('Converting PTKAirways result into image')
airways = lungs.BlankCopy();
airways = PTKGetImageFromAirwayResults(airway_results.AirwayTree, airways, false, ptk_main.Reporting);
airways.Title = [casename, ' / ', 'Segmentation', ' / ', 'PTKAirways'];

%% Process skeleton results into binary matrix as PTKimage
reporting.LogVerbose('Converting PTKAirwayCentreline result into image')
skeleton = lungs.BlankCopy();

new_image = zeros(skeleton.ImageSize, 'uint8');
new_image(skeleton.GlobalToLocalIndices(skeleton_results.CentrelinePoints)) = 1;
new_image(skeleton.GlobalToLocalIndices(skeleton_results.BifurcationPoints)) = 1;
skeleton.ChangeRawImage(new_image);

% sets origin of cropped volume to image volume
skeleton.SetVoxelToThis(skeleton_results.StartPoint, 1);
skeleton.Title = [casename, ' / ', 'Segmentation', ' / ', 'PTKAirwayCentreline'];

%% Export outputs compressed niftis
reporting.LogVerbose('Exporting all segmentation results as nifti files')
MimSaveAsNifti(image, savepath, [casename, '_PTKimage.nii'],reporting)
MimSaveAsNifti(lungs, savepath, [casename, '_PTKlungs.nii'], reporting)
MimSaveAsNifti(lobes, savepath, [casename, '_PTKlobes.nii'], reporting)
MimSaveAsNifti(airways, savepath, [casename, '_PTKairways.nii'], reporting)
MimSaveAsNifti(skeleton, savepath, [casename, '_PTKairwaycentreline.nii'],reporting)

%gzip niftifiles and delete originals
reporting.LogVerbose('Compressing output nifti files using gzip')
cd(savepath)
gzip('*_PTK*.nii')
% only deletes after originals gzipped
delete *_PTK*.nii
cd(initial_path)

%% save 2D previews of each segmentation
cd(savepath)
reporting.LogVerbose('Generate and save 2D quick previews of all segmentations')
% get CT image
lung_image = dataset.GetResult('PTKLungROI');

% lungs
image_viewer = PTKViewer(lung_image);
image_viewer.ViewerPanelHandle.OverlayImage = lungs;
image_viewer.ViewerPanelHandle.OverlayOpacity = 30;
PTKShow2DSlicesInOneFigure(image_viewer.ViewerPanelHandle,PTKImageOrientation.Axial, 20, reporting);
saveas(gcf,[casename, '_PTKlungs_quickresult'],'tiff')
close all

% airway s
image_viewer = PTKViewer(lung_image);
image_viewer.ViewerPanelHandle.OverlayImage = airways;
image_viewer.ViewerPanelHandle.OverlayOpacity = 30;
PTKShow2DSlicesInOneFigure(image_viewer.ViewerPanelHandle,PTKImageOrientation.Coronal, 60, reporting);
saveas(gcf,[casename, '_PTKairways_quickresult'],'tiff')
close all

% lobes
image_viewer = PTKViewer(lung_image);
image_viewer.ViewerPanelHandle.OverlayImage = lobes;
image_viewer.ViewerPanelHandle.OverlayOpacity = 30;
PTKShow2DSlicesInOneFigure(image_viewer.ViewerPanelHandle,PTKImageOrientation.Sagittal, 20, reporting);
saveas(gcf,[casename, '_PTKlobes_quickresult'],'tiff')
close all

cd(initial_path)

%% save snapshots of 3D models
cd(savepath)
reporting.LogVerbose('Generate and save 3D quick previews of all results')

% lungs
PTKVisualiseIn3D([], lungs, 4, false, reporting);
view(0,0)
camlight('right')
saveas(gcf,[casename, '_PTKlungs3D'],'png')
close all

% airways
PTKVisualiseIn3D([], airways, 0, true, reporting);
view(0,0)
camlight('right')
saveas(gcf,[casename, '_PTKairways3D'],'png')
close all

% skeleton
PTKVisualiseIn3D([], skeleton, 0, true, reporting);
view(0,0)
camlight('right')
saveas(gcf,[casename, '_PTKairwaycentreline3D'],'png')
close all

% lobes
PTKVisualiseIn3D([], lobes, 4, false, reporting);
view(0,0)
camlight('right')
saveas(gcf,[casename, '_PTKlobes3D'],'png')
close all

cd(initial_path)
reporting.Log('Complete')
