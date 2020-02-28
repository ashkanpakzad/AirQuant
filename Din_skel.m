% Skeletonisation algorithm implementation of D Jin et al. by Ashkan Pakzad
% on the 26th Feb 2020

%Getting the binary image
seg_name = 'github_demo_seg.nii.gz';
S = logical(niftiread(seg_name));

%% Distance Transform
St = ~S; % inverse
DTmap = bwdist(St, 'euclidean');

% Now get the local max in a 5x5 window everywhere.
localMaxImage = imdilate(DTmap, true(5)); % Whatever window size you want.

%% initialise


%% First iteration

