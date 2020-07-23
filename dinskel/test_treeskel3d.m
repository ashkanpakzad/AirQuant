%% testing treeskel 3d
prefix = '3D_sample';
suffix = '.tif';
n = 40;
object = zeros(512/4, 512/4, n);
% import noisy seg example
for i = 1:n
    j = i - 1;
    name = [prefix sprintf('%04d',j) suffix];
    t = imread(name);
    t = t(:,:,1); 
    t(t==255)=1;
    t = double(t);
    J = imresize(t, 0.25, 'nearest');
    object(:,:,i) = J;
end
object = logical(object);
% figure
% imagesc(proj)
% colormap gray
% hold on
inity=59; % swapped from plot
initx=63;
initz = 6;
init = [inity, initx, initz];
% plot(inity,initx,'.g')

%% run algo
tic;
Skel = TreeSkel3D(object, init, 3, 0.5, 0.5, 2, 0);
timeit = toc/60;

%% display against matlab skeletonise
figure
subplot(1,2,1)

skelpath = find(Skel);
[Y, X, Z] = ind2sub(size(object), skelpath);
patch(isosurface(object),'EdgeColor', 'none','FaceAlpha',0.3);
title('Implementation of Din et al. 2016')
hold on
plot3(X,Y,Z,'r.','MarkerSize',10)


subplot(1,2,2)
matlabskel = find(bwskel(logical(object)));
[Y, X, Z] = ind2sub(size(object), matlabskel);
patch(isosurface(object),'EdgeColor', 'none','FaceAlpha',0.3);
title('MATLAB''s implementation of Lee et al. 1994')
hold on
plot3(X,Y,Z,'r.','MarkerSize',10)

%%
se = strel('sphere',1);
boundaryimage = object - imerode(object, se);
figure
imshow(object(:,:,20), [])
figure
imshow(boundaryimage(:,:,35), [])