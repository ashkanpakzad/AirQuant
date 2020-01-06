%% Hyper Params
CT_name = 'github_demo_raw.nii.gz';
seg_name = 'github_demo_seg.nii.gz';
params = [];

%% Load data
meta = niftiinfo(CT_name);
CT = double(niftiread(meta));

%Getting the segmented image
S = logical(niftiread(seg_name));

%% create airway skeleton class
AS = AirwaySkel(CT, meta, S, params);

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

%% traverse one airway
AS = CreateAirwayImage(AS, 1);