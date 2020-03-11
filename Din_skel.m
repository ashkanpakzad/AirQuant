% Skeletonisation algorithm implementation of D Jin et al. by Ashkan Pakzad
% on the 26th Feb 2020

%Getting the binary image
seg_name = 'github_demo_seg.nii.gz';
S = logical(niftiread(seg_name));

%% create a projection to get a 2d binary image
proj = sum(S(:,:,86:90),3);
proj(proj~=0)=1;

imshow(proj,[])

%% Get Distance Transform
St = ~proj; % inverse
DTmap = bwdist(St, 'euclidean');
imshow(DTmap,[])

%% Create Central maximal disk look up table
T = table('Size',[0 4],'VariableTypes',{'double','double','double','double'});
T.Properties.VariableNames = {'radius', 'linear','X', 'Y'};
DTmapscan = DTmap;

candidatemax = 2; % to initialise
while candidatemax > 1
    % get current max
    [candidatemax, I] = max(DTmapscan,[],'all','linear');
    [Y, X] = ind2sub(size(DTmapscan), I);

    % check if inside an established CMB
    iscmb = 1;
    for i = 1:size(T)
        d = sqrt((X-T.X(i))^2 + (Y-T.Y(i))^2);
        if T.radius(i) > d || candidatemax == 0
            iscmb = 0;
            break % if this statement is true, then not a CMB
        end
    end
    
    % Add to cmb table
    if iscmb == 1
        newrow = {candidatemax, I, X, Y};
        T = [T; newrow];
    end
    DTmapscan(I) = 0;
    [candidatemax, I] = max(DTmapscan,[],'all','linear');
end

%% show all cmbs
figure
imagesc(DTmap)
colormap gray
hold on
for i = 1:size(T, 1)
circle(T.X(i), T.Y(i), T.radius(i));
end

CMBmap = zeros(size(DTmap));
CMBmap(T.linear) = 1;

figure
imagesc(CMBmap)
colormap gray

%% Computing LSF of all CMB
allLSF = zeros(size(T, 1),1);
for j = 1:size(T, 1)
py= T.X(j);
px= T.Y(j);

[Qx, Qy] = ndgrid(px-1:px+1, py-1:py+1);
nb = [Qx(:),Qy(:)];
nb(5,:) = [];

inequalityvals = zeros(length(neighbors),1);
for i = 1:length(neighbors)
    dist = abs(px-neighbors(i,1))+abs(py-neighbors(i,2));
    upper = DTmap(neighbors(i,1), neighbors(i,2)) - DTmap(px, py);
    lower = 0.5*(proj(px, py)+proj(neighbors(i,1), neighbors(i,2)))*dist;
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

%% Work out step cost
py= 285;
px= 248;
qy= 284;
qx= 248;

epsilon = 0.01;
LSF1 = sigfactor(px, py, proj, DTmap);
LSF2 = sigfactor(qx, qy, proj, DTmap);
dist = abs(px-qx)+abs(py-qy);

SC = dist/(epsilon+(LSF1+LSF2)^2);

%% find lowest step cost
% 
startx= 296;
starty= 273;
stopx=286;
stopy=197;

figure
h = imagesc(DTmap);
colormap gray
hold on

%plot start and stop
plot(startx,starty, 'g.')
plot(stopx,stopy, 'g.')

stop = 0;
j = 0;
totalSC = [];
while stop == 0
    if j == 0
        currentx = startx;
        currenty = starty;
    end
    j = j + 1;
    
    % plot
    plot(currentx, currenty, 'r.')
    
    % identify neighbors
    nb = neighbors(currentx, currenty);
    candidatecost = zeros(length(nb),1);
    for i = 1:length(nb)
        % compute step cost to each neighbor
        candidatecost(i) = stepcost(currentx,currenty,nb(i,1),nb(i,2),proj,DTmap);
    end
    % identify lowest cost neighbor
    [lowestSC,I] = min(candidatecost(:));
    
    % set new point as lowest cost neighbor.
    currentx = nb(I,1);
    currenty = nb(I,2);
    
    totalSC = [totalSC; lowestSC];
    
    if currentx == stopx && currenty == stopy
        stop = 1;
    end
    pause(0.5)
end



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
[Qx, Qy] = ndgrid(px-1:px+1, py-1:py+1);
neighbors = [Qx(:),Qy(:)];
neighbors(5,:) = [];

inequalityvals = zeros(length(neighbors),1);
for i = 1:length(neighbors)
    dist = abs(px-neighbors(i,1))+abs(py-neighbors(i,2));
    upper = DTmap(neighbors(i,1), neighbors(i,2)) - DTmap(px, py);
    lower = 0.5*(object(px, py)+object(neighbors(i,1), neighbors(i,2)))*dist;
    inequalityvals(i) = upper/lower;
end
output = abs(max(inequalityvals(:)));
if output > 0
    LSF = 1- output;
else
    LSF = 1;
end
end

function SC = stepcost(px,py,qx,qy,proj,DTmap)

epsilon = 0.01;
LSF1 = sigfactor(px, py, proj, DTmap);
LSF2 = sigfactor(qx, qy, proj, DTmap);
dist = abs(px-qx)+abs(py-qy);

SC = dist/(epsilon+(LSF1+LSF2)^2);
end

function nb = neighbors(px, py)
[Qx, Qy] = ndgrid(px-1:px+1, py-1:py+1);
nb = [Qx(:),Qy(:)];
nb(5,:) = [];
end