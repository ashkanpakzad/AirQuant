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

candidatemax = 1; % to initialise
while candidatemax > 0
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

%% functions
function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit, 'r');
hold off
end
