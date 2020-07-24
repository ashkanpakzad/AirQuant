function Skel = TreeSkel3D(object, init, thresh_min, thresh_multi, thresh_CMB, thresh_fill, debug)
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
Omarked = zeros(size(object));
Omarked(init(1), init(2), init(3)) = 1;
skelpath = []; % intialise skel path

%% Create Central maximal disk look up table
T = zeros(0,6);
%T VariableNames = {'radius', 'linear','X', 'Y', 'Z', 'LSF'};
DTmapscan = DTmap;

candidatemax = 1; % to initialise
counter = 1;
while candidatemax > 0
    counter = counter + 1;
    % get current max
    [candidatemax, I] = max(DTmapscan,[],'all','linear');
    [YI, XI, ZI] = ind2sub(size(DTmapscan), I);
    
    iscmb = 1;
    
    nb = neighbors_idx(XI, YI, ZI);
    nbind = sub2ind(size(DTmap),nb(2,:),nb(1,:),nb(3,:));
    neighbor_compare = 2*(DTmap(nbind)-candidatemax)./(object(I)+ object(nbind));
    
    if any(neighbor_compare >= 1)
        iscmb = 0;
        % if this statement is true, then not a CMB
    end
    
    % Add to cmb table
    if iscmb == 1
        newrow = [candidatemax, I, XI, YI, ZI, 0];
        T = [T; newrow];
    end
    DTmapscan(I) = 0;
    % used to identify if while condition met
    candidatemax = max(DTmapscan,[],'all','linear');
end


%%% Compute LSF of CMBs
allLSF = zeros(size(T, 1),1);
for j_lsf = 1:size(T, 1)
    py= T(j_lsf,3);
    px= T(j_lsf,4);
    pz = T(j_lsf,5);
    
    allLSF(j_lsf) = sigfactor(px, py, pz, object, DTmap);
end

T(:,6) = allLSF;

% identify strong CMBs
strongCMB = find(T(:,6) >= thresh_CMB);
%% Debugging CMB plots
if debug == 1
%     sparseCMB = strongCMB(1:100:end);
%     figure
%     patch(isosurface(object),'EdgeColor', 'none','FaceAlpha',0.3);
%     hold on
%     axis vis3d
%     for m = 1:size(sparseCMB, 1)
%         ball(T(sparseCMB(m),3), T(sparseCMB(m),4), ...
%             T(sparseCMB(m),5), T(sparseCMB(m),1));
%     end
    
    figure
    patch(isosurface(object),'EdgeColor', 'none','FaceAlpha',0.3);
    hold on
    plot3(T(:,3), T(:,4),T(:,5),'c.')
    plot3(T(strongCMB,3), T(strongCMB,4),T(strongCMB,5),'r.')
    legend('weak CMB', 'Strong CMB')
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

%%% while loop for as long as new branches found
branchit = 0;
while ~isequal(Omarked, object)
    branchit = branchit + 1;
    % detect potential subtree
    [subtrees, Tmarked] = subtree(object, Omarked, strongCMB);
    
    % update Omarked from potential subtree detector
    Omarked = (Omarked == 1 | Tmarked == 1);
    
    % checking here
%     figure
%     patch(isosurface(Omarked),'EdgeColor', 'none','FaceAlpha',0.3);
%     axis vis3d
%     axis([0 size(Omarked, 1) 0 size(Omarked, 2) 0 size(Omarked, 3)])
    
    if isempty(subtrees)
        break
    end
    
    % if potential branchs identified add branch
    for i_st = 1:length(subtrees)
        [furthestCMBX, furthestCMBY, furthestCMBZ] = furthestCMB(Omarked,T(strongCMB,2), subtrees{1,i_st});
        
        [pathX,pathY,pathZ,~] = mincostpath(...
            init(1), init(2), init(3), furthestCMBY, furthestCMBX,furthestCMBZ, object, g);
        % convert to linear indicies and extract branch
        newskelI = sub2ind(size(object), pathX, pathY, pathZ);
        Bi = setdiff(newskelI,skelpath);
        % check significance
        % add up LSF values along branch path to get LSFBi.
        Bi_LSF = zeros(size(Bi));
        for k = 1:length(Bi)
            [px, py, pz] = ind2sub(size(object), Bi(k));
            Bi_LSF(k) = sigfactor(px, py, pz, object, DTmap);
        end
        % compute DT value of CMB that branches the new skel to get DT_CMBv.
        if branchit ~= 1
            branch_v = newskelI(end-length(Bi),1);
            DT_branchv = DTmap(branch_v);
            % assume not sig
            Bi_sig = 0;
            % significant branch if this condition met
            sig_value = sum(Bi_LSF);
            adaptive_threshold = thresh_min+thresh_multi*(DT_branchv);
            if sig_value > adaptive_threshold
                Bi_sig = 1;
                %disp(int2str(sig_value))
            end
        else
            Bi_sig = 1;
        end
        % add branch to skel if sig
        if Bi_sig == 1
            skelpath = [skelpath; Bi];
        end
        % mark volume of new branch
        Bmarked = adaptivefill(Bi, object, DTmap);
        % update marked volume
        Omarked = (Omarked == 1 | Bmarked == 1);
%     
    end
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
        
        startind = sub2ind(size(object), sx, sy, sz);
        stopind = sub2ind(size(object), tx, ty, tz);
        
        % compute lowest cost path
        startnode = find(g.Nodes.PixelIndex == startind);
        stopnode = find(g.Nodes.PixelIndex == stopind);
        
        [P, cost] = shortestpath(g, startnode, stopnode);
        
        % get coords
        [pathX,pathY,pathZ] = ind2sub(size(object), g.Nodes.PixelIndex(P));
    end


    function Bmarked = adaptivefill(Bskel, object, DTmap)
        % Bskel = skeleton of branch only
        
        skel = zeros(size(object));
        skel(Bskel) = 1;
        % identify object but not skeleton voxels.
        nonskel = zeros(size(object));
        nonskel(object == 1 & skel == 0) = 1;
        
        % initialise with Dilation Scale map
        DSmap = zeros(size(object));
        DSmap(Bskel) = thresh_fill*DTmap(Bskel);
        DSmap(nonskel == 1) = -1*max(DSmap(:));
        
        % identify all p, part of object but not the skeleton
        allp = find(nonskel == 1);
        [allpx, allpy, allpz] = ind2sub(size(nonskel), allp);
        
        %copy DSmap
        DSmapprev = DSmap;
        DSmapnew = zeros(size(DSmap)); % holder!
        dsvals = zeros(nb_con,1);
        k_ds = 0;
        while ~isequal(DSmapnew, DSmapprev) % while there is no change
            % dummy on very first iteration
            if k_ds == 0
                DSmapnew = DSmap;
            end
            %update iteration
            k_ds = k_ds+1;
            % shift current to previous
            DSmapprev = DSmapnew;
            % loop through all voxels in object but not in skel
            for j_ds = 1:length(allp)
                
                px = allpx(j_ds);
                py = allpy(j_ds);
                pz = allpz(j_ds);
                
                % compare neighbors of voxel to current voxel.
                nb = neighbors_idx(px,py,pz);

                for i_ds = 1:nb_con
                    %dist = abs(px-nb(i_ds,1))+abs(py-nb(i_ds,2));
                    dist = sqrt((px-nb(1,i_ds))^2+(py-nb(2,i_ds))^2+(pz-nb(3,i_ds))^2);
                    dsvals(i_ds) = DSmapprev(nb(1,i_ds), nb(2,i_ds), nb(3,i_ds)) - dist;
                end

                DSmapnew(px,py,pz) = max(dsvals(:));
                
            end
        end
        
        Bmarked = (DSmapnew >= 0 & object == 1);
    end

    function [furthestCMBX, furthestCMBY, furthestCMBZ] = furthestCMB(Omarked, CMB_list, subtree_ind)
        % subtree_ind = cell2mat(CC_unmarked.PixelIdxList(B_potential(1)));
        % Omarked = Omarked;
        % CMB_list= T.linear;
        
        % subtree_ind = list of linear indicies in the subtree
        % Omarked = marked volume of skel algorithm
        % CMB_list = list linear indicies of Central Maximal Balls
        
        % get boundary voxels of Omarked
        se = strel('sphere',1);
        boundaryimage = Omarked - imerode(Omarked, se);
        bound_lin = find(boundaryimage == 1);
        [Yb, Xb, Zb] = ind2sub(size(boundaryimage), bound_lin);
        Omarked_bound = [Yb, Xb, Zb];
        
        % get subset of CMB in subtree
        commonCMBind = intersect(CMB_list, subtree_ind);
        
        % compute distances of all Omarked boundary to all CMB in T_i
        [commonCMBlinY, commonCMBlinX, commonCMBlinZ]= ind2sub(size(Omarked), commonCMBind);
        
        disti = zeros(length(commonCMBind),1);
        for i = 1:length(commonCMBind)
            distj = zeros(size(Omarked_bound,1),1);
            for j = 1:size(Omarked_bound,1)
                distj(j,1) = sqrt((commonCMBlinX(i,1) - Omarked_bound(j,1))^2 + (commonCMBlinY(i,1) - Omarked_bound(j,2))^2 + (commonCMBlinZ(i,1) - Omarked_bound(j,3))^2);
            end
            [disti(i,1)] = max(distj(:));
        end
        
        % extract max distance and return furthest CMB
        [~,furthestCMBcommon] = max(disti);
        furthestCMBX = commonCMBlinX(furthestCMBcommon);
        furthestCMBY = commonCMBlinY(furthestCMBcommon);
        furthestCMBZ = commonCMBlinZ(furthestCMBcommon);

    end

    function [subtrees, Tmarked] = subtree(object, Omarked, strongCMB)
        % subtract to get unmarked O
        Tmarked = Omarked;
        Ounmarked = zeros(size(Tmarked));
        Ounmarked(object == 1 & Tmarked == 0) = 1;
        
        % identify potential branches by connectivity
        CC_unmarked = bwconncomp(Ounmarked);
        
        % identify if any branches are significant
        % loop through all unmarked components
        B_potential = [];
        for i = 1:length(CC_unmarked.PixelIdxList)
            CMBinCC = 0;
            % identify if any strong CMB in unmarked components
            for j = 1:length(strongCMB)
                if ~isempty(find(cell2mat(CC_unmarked.PixelIdxList(i))==T(strongCMB(j),2),1))
                    % mark CMB in current CC
                    CMBinCC = 1;
                end
            end
            if CMBinCC == 1
                B_potential = [B_potential; i];
            else
                % if no strong CMB mark components
                Tmarked(cell2mat(CC_unmarked.PixelIdxList(i))) = 1;
            end
        end
        subtrees = CC_unmarked.PixelIdxList(1, B_potential);
    end

    function h = ball(x,y,z,r)
        hold on
        [X,Y,Z] = sphere;
        X2 = X * r + x;
        Y2 = Y * r + y;
        Z2 = Z * r + z;
        surf(X2,Y2,Z2)
    end
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