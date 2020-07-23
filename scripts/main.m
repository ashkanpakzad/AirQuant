%% Hyper Params
CT_name = 'github_demo_raw.nii.gz';
seg_name = 'github_demo_seg.nii.gz';
params = [];

% Load data
meta = niftiinfo(CT_name);
CT = double(niftiread(meta));

%Getting the segmented image
S = logical(niftiread(seg_name));

% create airway skeleton class
AS = AirwaySkel(CT, meta, S, params);

%% plot
figure
PlotTree(AS)
hold on
D = logical(niftiread('github_demo_dital_point.nii'));
se = strel('sphere',4);
D = imdilate(D, se);
patch(isosurface(D),'FaceColor','none','EdgeColor',[0 1 0]);


%% traverse all airway
tic;
AS = AirwayImageAll(AS)
time = toc;
save('CompletedAirway2')

%% compute area of all airways
AS = FindFWHMall(AS)

%% Construct taper rate path to compute
[logtaperrate, cum_arclength, cum_area, path] = ConstructTaperPath(obj, terminal_link_idx); 

%% Construct spline, find traversed image for first point
spline = ComputeSpline(AS, 1);
spline_para_limit = spline.breaks(end);
spline_points = 0:AS.spline_sampling_interval:spline_para_limit;

% loop along spline
TransAirwayImage = zeros(133,133,length(spline_points));
for i = 1:4
    % * Compute Normal Vector per spline point
    [normal, CT_point] = AirwaySkel.ComputeNormal(spline, spline_points(i));
    % * Interpolate Perpendicular Slice per spline point
    TransAirwayImage(:,:,i) = InterpolateCT(AS, normal, CT_point);
end

%% have a go at slider
sliderchanging(AS, 1)
%%
% function sliderchanging(AS, link)
% % Create figure window and components
% 
% figure
% sld = uicontrol('style','slider','units','pixel','position',[20 20 300 20]);
% 
%                'ValueChangingFcn',@(sld,event) sliderMoving(event,AS,link));
% sld.Limits = [1 size(AS.TraversedImage{link,1}, 3)];
% sld.Value = 1;
% 
% end
% 
% % Create ValueChangingFcn callback
% function sliderMoving(event,AS,link)
% val = round(event.Value);
% imshow(AS.TraversedImage{link,1}(:,:,val), [])
% title(['Traversed link = ', num2str(link),...
%     ' arclength = ', num2str(AS.spline_sampling_interval*val - 0.25)])
% drawnow;
% end
% 
% function myslider
% x = 1:10;
% hplot = plot(x,0*x);
% 
% 
% 
% 


