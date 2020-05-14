function Skel = TreeSkel_terminal(object, terminal, debug)
% Object = binary image
% init = point to start skeletonisation from.

% default threshmin = 3
% default threshmulti = 0.5
% default thresh_CMB = 0.5
% default thresh_fill = 1.5

% TODO: catch error where init is outside of object...

% set up record of neighbor positions
nb_con = 26; % pixel connectivity

inequalityvals = zeros(nb_con,1); % reserve memory for LSF computing.

%% get DT map & Omarked
DTmap = bwdist(~object, 'euclidean');
skelpath = []; % intialise skel path
init = zeros(1,3);
% points marked by label of 2 is the trachea terminus.
[init(2), init(1), init(3)] = ind2sub(size(object), find(terminal == 2));
%terminal(init(1), init(2), init(3)) = 0; % make sure it doesnt appear in T.

%% Create list of terminal points
T = find(terminal == 1);
% identify where trachea terminal is and swap it to 1st position
% init_lin = sub2ind(size(object), init(2), init(1), init(3));
% tmp_init_idx = (T == init_lin); % Trachea point
% tmp_1 = T(1);
% T(1) = init_lin;
% T(tmp_init_idx) = tmp_1;

% get sub indices of all terminals and add columns.
[T(:,3), T(:,2), T(:,4)] = ind2sub(size(object), T);

%T VariableNames = {'radius', 'linear','X', 'Y', 'Z', 'LSF'};

%% Debugging CMB plots
if debug == 1
    figure
    patch(isosurface(object),'EdgeColor', 'none','FaceAlpha',0.3);
    hold on
    plot3(init(1), init(2), init(3),'r.')
    plot3(T(:,2), T(:,3),T(:,4),'c.')
    legend('object','Trachea Terminus', 'Airway Terminus')
    axis vis3d
    error('Debug mode called.') % break out of function
end

%% Compute loss graph
g = binaryImageGraph3(object);

% extract graph information into matrices for performance.
p_all = g.Edges.EndNodes(:,1);
q_all = g.Edges.EndNodes(:,2);
nodex = g.Nodes.x;
nodey = g.Nodes.y;
nodez = g.Nodes.z;
% compute step cost across grid

weights = zeros(size(p_all));
for iedge = 1:height(g.Edges)
    pnode = p_all(iedge);
    qnode = q_all(iedge);
    py = nodex(pnode);
    px = nodey(pnode);
    pz = nodez(pnode);
    qy = nodex(qnode);
    qx = nodey(qnode);
    qz = nodez(qnode);
    weights(iedge) = stepcost(px,py,pz,qx,qy,qz,object,DTmap);
end

g.Edges.Weight = weights;

%% while loop for as long as terminals not all added
branchit = 0;
for b = 1:size(T,1)
    branchit = branchit + 1;
    % identify furthest terminal

    [pathX,pathY,pathZ,~] = mincostpath(...
        init(1), init(2), init(3), T(b,3), T(b,2),T(b,4), object, g);
    % convert to linear indicies and extract branch
    newskelI = sub2ind(size(object), pathX, pathY, pathZ);
    Bi = setdiff(newskelI,skelpath);
    
    % compute DT value of CMB that branches the new skel to get DT_CMBv.
    skelpath = [skelpath; Bi];
end

Skel = zeros(size(object));
Skel(skelpath) = 1;

    function LSF = sigfactor(px, py, pz, object, DTmap)
        nb = neighbors_idx(px,py,pz);
        for i = 1:nb_con
            dist = abs(px-nb(1,i))+abs(py-nb(2,i))+abs(pz-nb(3,i));
            upper = DTmap(nb(1,i), nb(2,i), nb(3,i)) - DTmap(px, py ,pz);
            lower = 0.5*(object(px,py,pz)+object(nb(1,i), nb(2,i), nb(3,i)))*dist;
            inequalityvals(i) = upper/lower;
        end
        output = abs(max(inequalityvals(:)));
        if output > 0
            LSF = 1 - output;
        else
            LSF = 1;
        end
    end

    function SC = stepcost(px,py,pz,qx,qy,qz,object,DTmap)
        
        epsilon = 0.01;
        LSF1 = sigfactor(px, py, pz, object, DTmap);
        LSF2 = sigfactor(qx, qy, qz, object, DTmap);
        dist_sc = abs(px-qx)+abs(py-qy)+abs(pz-qz);
        
        SC = dist_sc/(epsilon+(LSF1+LSF2)^2);
    end

    function [pathX,pathY,pathZ,cost] = mincostpath(sx, sy, sz, tx, ty,tz, object, g)
        
        startind = sub2ind(size(object), sy, sx, sz);
        stopind = sub2ind(size(object), tx, ty, tz);
        
        % compute lowest cost path
        startnode = find(g.Nodes.PixelIndex == startind);
        stopnode = find(g.Nodes.PixelIndex == stopind);
        
        [P, cost] = shortestpath(g, startnode, stopnode);
        
        % get coords
        [pathX,pathY,pathZ] = ind2sub(size(object), g.Nodes.PixelIndex(P));
    end

%     function [terminus] = furthestT()
%         
%         %* compute distances of skel points to terminals
%         
%         % Convert T linear indices to sub.
%         [commonCMBlinY, commonCMBlinX, commonCMBlinZ]= ind2sub(size(Omarked), commonCMBind);
%         
%         % Convert skeleton linear indices to sub.
%         
%         % Get distances between all points.
%         
%         disti = zeros(length(commonCMBind),1);
%         for i = 1:length(commonCMBind)
%             distj = zeros(size(Omarked_bound,1),1);
%             for j = 1:size(Omarked_bound,1)
%                 distj(j,1) = sqrt((commonCMBlinX(i,1) - Omarked_bound(j,1))^2 + (commonCMBlinY(i,1) - Omarked_bound(j,2))^2 + (commonCMBlinZ(i,1) - Omarked_bound(j,3))^2);
%             end
%             [disti(i,1)] = max(distj(:));
%         end
%         
%         % extract max distance and return furthest CMB
%         [~,furthestCMBcommon] = max(disti);
%         furthestCMBX = commonCMBlinX(furthestCMBcommon);
%         furthestCMBY = commonCMBlinY(furthestCMBcommon);
%         furthestCMBZ = commonCMBlinZ(furthestCMBcommon);
%     end

    function nb = neighbors_idx(i, j, k)
        nb = [i-1 j-1 k-1; ...
            i-1 j k-1; i-1 j+1 k-1; i j-1 k-1; ...
            i j k-1; i j+1 k-1; i+1 j-1 k-1; ...
            i+1 j k-1; i+1 j+1 k-1; i-1 j-1 k; ...
            i-1 j k; i-1 j+1 k; i j-1 k; ...
            i j+1 k; i+1 j-1 k; ...
            i+1 j k; i+1 j+1 k; i-1 j-1 k+1; ...
            i-1 j k+1; i-1 j+1 k+1; i j-1 k+1; ...
            i j k+1; i j+1 k+1; i+1 j-1 k+1; ...
            i+1 j k+1; i+1 j+1 k+1;]';
    end
end