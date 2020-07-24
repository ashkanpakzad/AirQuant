
%% Testing on real segmentation...
addpath(genpath('/Users/apakz/PhD/normal_data/'))
casename = 'N3';
seg_name = [casename, '_lumen.nii.gz'];
skel_name = [casename, '_skel.nii.gz'];
PTKskel = logical(niftiread(skel_name));

terminals_name = [casename, '_terminal.nii.gz'];
%Getting the segmented image and manual skeleton
S = logical(niftiread(seg_name));
terms = niftiread(terminals_name);
%%
tic;
Skel = TreeSkel_terminal(S, terms, 0);
timeit = toc/60;
%% display against matlab skeletonise
figure
patch(isosurface(S),'EdgeColor', 'none','FaceAlpha',0.3);
title('N3 - Implementation of Jin et al. 2016 with PTK Terminus points')
hold on
isosurface(Skel)
axis vis3d
view(80,0)

%%
figure
Leeskel = bwskel(logical(S));
patch(isosurface(S),'EdgeColor', 'none','FaceAlpha',0.3);
title('N3 - MATLAB''s implementation of Lee et al. 1994')
hold on
isosurface(Leeskel)
axis vis3d
view(80,0)

figure
patch(isosurface(S),'EdgeColor', 'none','FaceAlpha',0.3);
title('N3 - PTK''s implementation of Palagyi 2006')
hold on
isosurface(PTKskel)
axis vis3d
view(80,0)

%% save skeleton

%%
se = strel('sphere',1);
boundaryimage = object - imerode(object, se);
figure
imshow(object(:,:,20), [])
figure
imshow(boundaryimage(:,:,35), [])