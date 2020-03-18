% Skeletonisation algorithm implementation of D Jin et al. by Ashkan Pakzad
% on the 26th Feb 2020

%Getting the binary image
seg_name = 'github_demo_seg.nii.gz';
S = logical(niftiread(seg_name));

%% create a projection to get a 2d binary image
proj = sum(S(:,:,86:90),3);
proj(proj~=0)=1;
proj = proj(150:300,220:330);
imshow(proj,[])

%% Get Distance Transform
St = ~proj; % inverse
DTmap = bwdist(St, 'euclidean');
imshow(DTmap,[])

%% Create Central maximal disk look up table
T = table('Size',[0 4],'VariableTypes',{'double','double','double','double'});
T.Properties.VariableNames = {'radius', 'linear','X', 'Y'};
DTmapscan = DTmap;

candidatemax = 1; % to initialise

B = bwboundaries(proj); % find boundaries
B = B{1,1};
while candidatemax > 0
    % get current max
    [candidatemax, I] = max(DTmapscan,[],'all','linear');
    [Y, X] = ind2sub(size(DTmapscan), I);

    iscmb = 1;
    
    % check if neighbor has a larger value
    nb = neighbors(X,Y);
    nbind = sub2ind(size(DTmap),nb(:,2),nb(:,1));
    nbDT = DTmap(nbind);
    neighbor_compare = 2*(nbDT-candidatemax)./(proj(I)+ proj(nbind));
    if any(neighbor_compare >= 1)
        iscmb = 0;
        % if this statement is true, then not a CMB
    end
    
    % check if inside an established CMB
%     if iscmb == 1
%         for i = 1:size(T)
%             d = sqrt((X-T.X(i))^2 + (Y-T.Y(i))^2);
%             if candidatemax == 0 || T.radius(i) > d + candidatemax
%                 iscmb = 0;
%                 break % if this statement is true, then not a CMB
%             end
%         end
%     end
    
    % Add to cmb table
    if iscmb == 1
        newrow = {candidatemax, I, X, Y};
        T = [T; newrow];
    end
    DTmapscan(I) = 0;
    [candidatemax, I] = max(DTmapscan,[],'all','linear');
end

%% Computing LSF of all CMB
allLSF = zeros(size(T, 1),1);
for j = 1:size(T, 1)
    py= T.X(j);
    px= T.Y(j);
    
    [Qx, Qy] = ndgrid(px-1:px+1, py-1:py+1);
    nb = [Qx(:),Qy(:)];
    nb(5,:) = [];
    
    inequalityvals = zeros(length(nb),1);
    for i = 1:length(nb)
        dist = abs(px-nb(i,1))+abs(py-nb(i,2));
        upper = DTmap(nb(i,1), nb(i,2)) - DTmap(px, py);
        lower = 0.5*(proj(px, py)+proj(nb(i,1), nb(i,2)))*dist;
        inequalityvals(i) = upper/lower;
    end
    output = abs(max(inequalityvals(:)));
    if output > 0
        LSF = 1- output;
    else
        LSF = 1;
    end
    allLSF(j) = LSF;
end
T.LSF = allLSF;

%% show all cmbs
figure
subplot(1,2,1)
imagesc(proj)
colormap gray
hold on
for i = 1:size(T, 1)
circle(T.X(i), T.Y(i), T.radius(i));
end

subplot(1,2,2)
imagesc(DTmap)
colormap gray
hold on
plot(B(:,2), B(:,1), 'g', 'LineWidth', 2)
plot(T.X, T.Y,'c.')
strongCMB = find(T.LSF >= 0.5);
plot(T.X(strongCMB), T.Y(strongCMB),'r.')
legend('Object Boundary', 'weak CMB', 'Strong CMB')


%% Work out step cost
px= 50;
py= 50;
qx= 62;
qy= 64;

epsilon = 0.01;
LSF1 = sigfactor(px, py, proj, DTmap);
LSF2 = sigfactor(qx, qy, proj, DTmap);
dist = abs(px-qx)+abs(py-qy);

SC = dist/(epsilon+(LSF1+LSF2)^2);

%% %% find lowest step cost by computing all steps
% 
startx= 72;
starty= 75;
stopx=67;
stopy=50;

startind = sub2ind(size(proj), starty, startx);
stopind = sub2ind(size(proj), stopy, stopx);

% convert image grid to graph
g = binaryImageGraph(proj);
%plotImageGraph(g)

% compute step cost across grid
for i = 1:height(g.Edges)
    pnode = g.Edges.EndNodes(i,1);
    qnode = g.Edges.EndNodes(i,2);
    px = g.Nodes.x(pnode);
    py = g.Nodes.y(pnode);
    qx = g.Nodes.x(qnode);
    qy = g.Nodes.y(qnode);
    g.Edges.Weight(i) = stepcost(px,py,qx,qy,proj,DTmap);
end

% compute lowest cost path
startnode = find(g.Nodes.PixelIndex == startind);
stopnode = find(g.Nodes.PixelIndex == stopind);

[P, cost] = shortestpath(g, startnode, stopnode);

% display
% get image of path
path = zeros(size(proj));
path(g.Nodes.PixelIndex(P)) = 1;

% get coords
[pathX,pathY] = ind2sub(size(proj), g.Nodes.PixelIndex(P));

%% compute mincostpath by function
sx= 121;
sy= 71;
tx=37;
ty=31;
[pathX,pathY,cost] = mincostpath(sx, sy, tx, ty, proj, DTmap);

figure
imagesc(DTmap)
colormap gray
hold on
plot(pathY,pathX, 'y.')
%plot start and stop
plot(sy,sx, 'g.')
plot(ty,tx, 'g.')

%% intialise and find furthest CMB 

sx= 121;
sy= 71;

% //identify furthest CMB
Sxy = ones(height(T),2);
Sxy(:,1) = sx;
Sxy(:,2) = sy;

% compute euclidean distances
CMB_dist = sqrt((Sxy(:,2)-T.Y).^2 + (Sxy(:,1)-T.X).^2);
[~, CMB_far_I] = max(CMB_dist);
ty=T.X(CMB_far_I);
tx=T.Y(CMB_far_I);

% // compute min cost path
[pathX,pathY,cost] = mincostpath(sx, sy, tx, ty, proj, DTmap);

% // display result
figure
imagesc(DTmap)
colormap gray
hold on
plot(pathY,pathX, 'y.')
%plot start and stop
plot(sy,sx, 'g.')
plot(ty,tx, 'g.')

%% add branch to skel
skel = zeros(size(proj));

pathI = sub2ind(size(skel), pathX, pathY);
skel(pathI) = 1;

figure
imagesc(skel)
colormap gray

%% dialate fill
% identify object but not skeleton voxels.
nonskel = zeros(size(proj));
nonskel(proj == 1 & skel == 0) = 1;

% initialise with Dilation Scale map
DSmap = zeros(size(proj));
DSmap(pathI) = DTmap(pathI);
DSmap(nonskel == 1) = -1*max(DSmap(:));

% identify all p part of object but not the skeleton
allp = find(nonskel == 1);
[allpx, allpy] = ind2sub(size(nonskel), allp);

%copy DSmap
DSmapprev = DSmap;
DSmapnew = zeros(size(DSmap)); % holder!
k = 0;
while ~isequal(DSmapnew, DSmapprev) % while there is no change
    % dummy on very first iteration
    if k == 0
        DSmapnew = DSmap;
    end
    %update iteration
    k = k+1;
    % shift current to previous
    DSmapprev = DSmapnew;
    % loop through all voxels in object but not in skel
    for j = 1:length(allp)
        
        px = allpx(j);
        py = allpy(j);
        
        % identify neighbor
        nb = neighbors(px, py);
        
        dsvals = zeros(length(nb),1);
        for i = 1:length(nb)
            %dist = abs(px-nb(i,1))+abs(py-nb(i,2));
            dist = sqrt((px-nb(i,1))^2+(py-nb(i,2))^2);
            dsvals(i) = DSmapprev(nb(i,1), nb(i,2)) - dist;
        end
        % update DSmap
        DSmapnew(px,py) = max(dsvals(:));
    end
end

dilatefill = (DSmapnew >= 0 & proj == 1);
figure
imagesc(dilatefill)
colormap gray
hold on
plot(B(:,2), B(:,1), 'g', 'LineWidth', 2)

%% functions
function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit, 'r');
hold off
end

function LSF = sigfactor(px, py, object, DTmap)
nb = neighbors(px, py);

inequalityvals = zeros(length(nb),1);
for i = 1:length(nb)
    dist = abs(px-nb(i,1))+abs(py-nb(i,2));
    DTval = DTmap(px, py);
    upper = DTmap(nb(i,1), nb(i,2)) - DTval;
    lower = 0.5*(object(px, py)+object(nb(i,1), nb(i,2)))*dist;
    inequalityvals(i) = upper/lower;
end
output = abs(max(inequalityvals(:)));
if output > 0
    LSF = 1- output;
else
    LSF = 1;
end
end

function SC = stepcost(px,py,qx,qy,object,DTmap)

epsilon = 0.01;
LSF1 = sigfactor(px, py, object, DTmap);
LSF2 = sigfactor(qx, qy, object, DTmap);
dist = abs(px-qx)+abs(py-qy);

SC = dist/(epsilon+(LSF1+LSF2)^2);
end

function nb = neighbors(px, py)
[Qx, Qy] = ndgrid(px-1:px+1, py-1:py+1);
nb = [Qx(:),Qy(:)];
nb(5,:) = [];
end

function [pathX,pathY,cost] = mincostpath(sx, sy, tx, ty, object, DTmap)

startind = sub2ind(size(object), sx, sy);
stopind = sub2ind(size(object), tx, ty);

% convert image grid to graph
g = binaryImageGraph(object);
%plotImageGraph(g)

% compute step cost across grid
for i = 1:height(g.Edges)
    pnode = g.Edges.EndNodes(i,1);
    qnode = g.Edges.EndNodes(i,2);
    py = g.Nodes.x(pnode);
    px = g.Nodes.y(pnode);
    qy = g.Nodes.x(qnode);
    qx = g.Nodes.y(qnode);
    g.Edges.Weight(i) = stepcost(px,py,qx,qy,object,DTmap);
end

% compute lowest cost path
startnode = find(g.Nodes.PixelIndex == startind);
stopnode = find(g.Nodes.PixelIndex == stopind);

[P, cost] = shortestpath(g, startnode, stopnode);

% get coords
[pathX,pathY] = ind2sub(size(object), g.Nodes.PixelIndex(P));
end