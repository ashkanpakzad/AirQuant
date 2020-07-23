% test treeskel
% load in binary object
seg_name = 'github_demo_seg.nii.gz';
S = logical(niftiread(seg_name));

%create a projection to get a 2d binary image
proj = sum(S(:,:,86:90),3);
proj(proj~=0)=1;
proj = proj(150:300,220:330);

%%
inity=54;
initx=28;

Skel = TreeSkel(proj, initx, inity);

%% import noisy seg example
t = imread('noisyairwayseg.tiff');
t = t(:,:,1);
t(t==255)=1;
t = double(t);
J = imresize(t, 0.5, 'nearest');
proj = logical(J);
figure
imagesc(proj)
colormap gray
hold on
inity=116; % swapped from plot
initx=34;
plot(inity,initx,'.g')

%% run algo
tic;
Skel = TreeSkel(proj, initx, inity);
timeit = toc/60;

%% display against matlab skeletonise
figure

subplot(1,2,1)
skelpath = find(Skel);
[Y, X] = ind2sub(size(proj), skelpath);
imagesc(proj)
title('Implementation of Din et al. 2016')
colormap gray
hold on
plot(X,Y,'r.','MarkerSize',10)

subplot(1,2,2)
matlabskel = find(bwskel(logical(proj)));
[Y, X] = ind2sub(size(proj), matlabskel);
imagesc(proj)
title('MATLAB''s implementation of Lee et al. 1994')
colormap gray
hold on
plot(X,Y,'r.','MarkerSize',10)
