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

