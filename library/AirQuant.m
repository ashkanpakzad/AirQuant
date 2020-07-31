classdef AirQuant
    properties
        CT % CT image
        CTinfo % CT metadata
        seg % binary airway segmentation
        % CT Properties/resampling params
        physical_plane_length = 40;% check methods
        physical_sampling_interval = 0.3;% check methods
        spline_sampling_interval = 0.25;% check methods
        % Ray params
        num_rays = 50; % check methods
        ray_interval = 0.2; % check methods
        skel % skeleton based on segementation
        savename % filename to load/save
    end
    properties (SetAccess = private)
        % Graph Properties
        Gadj % undirected Graph Adjacency matrix
        Gnode % Graph node info
        Glink % Graph edge info
        Gdigraph % digraph object
        trachea_node % node corresponding to top of trachea
        trachea_path % edges that form a connect subgraph above the carina
        carina_node % node that corresponds to the carina
        % Resampled image slices along graph paths/airway segments.
        TraversedImage % perpendicular CT slices of all airways
        TraversedSeg % perpendicular segmentation slices of all airways
        arclength % corresponding arclength measurement traversed slices
        FWHMesl % FWHMesl algorithm results for every airway
        Specs % Airway specs
    end
    
    methods
        %% INITIALISATION METHODS
        function obj = AirQuant(CTimage, CTinfo, segimage, skel, savename, params)
            % Initialise the AirQuant class object.
            % if using default settings, do not provide params argument.
            
            % Can provide just file name, if analyses previously complete.
            if nargin == 1
                savename = CTimage;
                obj.savename = savename;
            end
            
            if isfile(savename)
                load(savename)
                obj.savename = savename;
            else
                obj.CT = CTimage;
                % TODO: consider preprocess segmentation to keep largest
                % connected component.
                % ensure no holes in segmentation
                obj.seg = imfill(segimage,'holes');
                obj.CTinfo = CTinfo;
                % set params
                if nargin > 5
                    obj.physical_plane_length = params.physical_plane_length;
                    obj.physical_sampling_interval = params.physical_sampling_interval;
                    obj.spline_sampling_interval = params.spline_sampling_interval;
                    obj.num_rays = params.num_rays;
                    obj.ray_interval = params.ray_interval;
                end
                % graph airway skeleton
                obj = GenerateSkel(obj,skel);
                % Identify trachea
                obj = FindTrachea(obj);
                % Convert into digraph
                obj = AirwayDigraph(obj);
                % Identify Carina
                obj = FindCarina(obj);
                % Identify paths that belong to trachea
                obj = FindTracheaPaths(obj);
                % classify airway generations
                obj = ComputeAirwayGen(obj);
                % set up empty cell for traversed CT and segmentation
                obj.TraversedImage = cell(length(obj.Glink),1);
                obj.TraversedSeg = cell(length(obj.Glink),1);
                % set up empty specs doubles
                obj.arclength = cell(length(obj.Glink),1);
                obj.Specs = struct();
                % set up empty cell for recording raycast/fwhmesl method
                obj.FWHMesl = cell(length(obj.Glink),3);
                % save class
                obj.savename = savename;
                save(obj)
            end
        end
        
        
        function obj = GenerateSkel(obj, skel)
            % Generate the airway skeleton
            % if skeleton not provided, generate one.
            if isempty(skel)
                obj.skel = Skeleton3D(obj.seg);
            else
                obj.skel = skel;
            end
            % create graph from skeleton.
            [obj.Gadj,obj.Gnode,obj.Glink] = Skel2Graph3D(obj.skel,0);
            
        end
        
        
        function obj = FindTrachea(obj)
            % Identify the Trachea node. FindFWHMall
            % Assumes trachea fully segmented and towards greater Z.
            % Smoothen
            [~, obj.trachea_node] = max([obj.Gnode.comz]);
            %obj.trachea_path = obj.Gnode(obj.trachea_node).links;
        end
        
        
        function obj = AirwayDigraph(obj)
            % Converts the output from skel2graph into a digraph tree
            % network, such that there are no loops and edges point from
            % the trachea distally outwards.
            
            % Create digraph with edges in both directions, loop through
            % and remove if found later in the BF search.
            G = digraph(obj.Gadj);
            % BF search from carina node to distal.
            node_discovery = bfsearch(G,obj.trachea_node);
            % half of the edges will be removed
            removal = zeros(height(G.Edges)/2 , 1);
            j = 1;
            for i = 1:height(G.Edges)
                % identify the rank in order of discovery
                rank1 = find(~(node_discovery-G.Edges.EndNodes(i,1)));
                rank2 = find(~(node_discovery-G.Edges.EndNodes(i,2)));
                if rank1 > rank2
                    % record if not central-->distal
                    removal(j) = i;
                    j = j + 1;
                end
            end
            removal(removal == 0) = [];
            G = rmedge(G,removal);
            
            % Need to ensure directions in Glink are also in the correct
            % direction. Swap over if not.
            for i = 1:length(obj.Glink)
                glink_nodes = [obj.Glink(i).n1, obj.Glink(i).n2];
                % if not central-->distal
                if ~ismember(glink_nodes, G.Edges.EndNodes, 'rows')
                    % Swap round node assignments
                    obj.Glink(i).n1 = glink_nodes(2);
                    obj.Glink(i).n2 = glink_nodes(1);
                    % Swap round path order
                    obj.Glink(i).point = fliplr(obj.Glink(i).point);
                end
            end
            
            % Copy over Glink to create digraph with same edge order
            edges = [[obj.Glink.n1]', [obj.Glink.n2]'];
            weights = zeros(length(obj.Glink),1);
            for i = 1:length(obj.Glink)
                weights(i) = length(obj.Glink(i).point);
            end
            labels = [1:length(obj.Glink)]';
            Edgetable = table(edges,weights,labels,'VariableNames',{'EndNodes', 'Weight', 'Label'});
            
            obj.Gdigraph = digraph(Edgetable);
            % add node properties from Gnode
            obj.Gdigraph.Nodes.comx(:) = [obj.Gnode(:).comx];
            obj.Gdigraph.Nodes.comy(:) = [obj.Gnode(:).comy];
            obj.Gdigraph.Nodes.comz(:) = [obj.Gnode(:).comz];
            obj.Gdigraph.Nodes.ep(:) = [obj.Gnode(:).ep];
            obj.Gdigraph.Nodes.label(:) = [1:length(obj.Gnode)]';
            
        end
        
        function obj = FindCarina(obj)
            % Finds the carina node by analysis of the directed graph
            % produced by a breadth first search from the trachea node.
            [~, obj.carina_node] = max(centrality(obj.digraph,'outcloseness'));
        end
        
        function obj = FindTracheaPaths(obj)
            % Check if there are any edges between carina and trachea
            % due to skeletonisation error (trachea is prone to
            % skeletonisation errors)
            
            % identify outgoing nodes from carina
            [eid, nid] = outedges(obj.Gdigraph,obj.carina_node);
            % identify their importance and sort in that order
            outcloseness = centrality(obj.digraph,'outcloseness');
            [~, I_sort] = sort(outcloseness(nid),'descend');
            % remove most important === bronchi.
            G = rmedge(obj.Gdigraph,eid(I_sort(1:2)));
            % identify smallest connected graph corresponding to trachea.
            [bin,binsize] = conncomp(G,'Type','weak');
            idx = binsize(bin) == min(binsize);
            % idx corresponds to nodes in trachea group
            SG = subgraph(G, idx);
            % identify labels
            obj.trachea_path=table2array(SG.Edges(:, {'Label'}));
        end
        
        function obj = ComputeAirwayGen(obj)
            % loop through graph nodes
            G = obj.Gdigraph;
            gens = zeros(length(obj.Glink),1);
            gens(obj.trachea_path) = 0;
            for i = 1:height(G.Nodes)
                if i == obj.trachea_path
                    break
                end
                % get path from carina to current node
                path = shortestpath(G, obj.carina_node, i);
                % edges out of carina, gen = 1
                generation = length(path);
                edge_idx = outedges(G, i);
                % assign generation from number of nodes traversed to out
                % edges.
                gens(G.Edges.Label(edge_idx)) = generation;
            end
            % add gens field to glink
            gens = num2cell(gens);
            [obj.Glink(:).generation] = gens{:};
        end
        
        function obj = ComputeAirwayLobes(obj)
            % Get airway graph
            G = obj.Gdigraph;
            % Set up lobe store vector
            lobes = cell(length(obj.Glink),1);
            % Assign trachea to itself
            lobes{obj.trachea_path} = 'B';
            
            % % Identify left and right lung nodes first
            lungN = successors(G, obj.carina_node);
            if obj.Gnode(lungN(1)).comx > obj.Gnode(lungN(2)).comx
                leftN = lungN(1);
                rightN = lungN(2);
            else 
                leftN = lungN(2);
                rightN = lungN(1);
            end
            
            % assign labels to major bronchi
            classedgenonlobes(obj.carina_node, leftN, 'B')
            classedgenonlobes(obj.carina_node, rightN, 'B')
                        
            % % Identify node of upper-lingular and lower left lobe 'LL'
            LlungN = successors(G, leftN);
            
            if obj.Gnode(LlungN(1)).comz > obj.Gnode(LlungN(2)).comz
                LUL_L = LlungN(1);
                LLLN = LlungN(2);
            else
                LUL_L = LlungN(2);
                LLLN = LlungN(1);
            end
            
            % Assign branches of the lower left lobe
            classedgelobes(LLLN, 'LL')
            classedgelobes(LUL_L, 'LU')

            % assign branch from left lobe divider to upper lobe node
            classedgenonlobes(leftN, LUL_L, 'LU')
            
            % identify upper lobe and lingular
            LUL_LN = successors(G, LUL_L);
            
            if obj.Gnode(LUL_LN(1)).comz > obj.Gnode(LUL_LN(2)).comz
                LULN = LUL_LN(1);
                LN = LUL_LN(2);
            else
                LULN = LUL_LN(2);
                LN = LUL_LN(1);
            end
            
            % assign branches of the upper left lobe.
            classedgelobes(LULN, 'LU')
            classedgelobes(LN, 'lin')
                        
            % % Identify right upper lobe
            RlungN = successors(G, rightN);
            [~, I] = max([obj.Gnode(RlungN).comz]);
            RULN = RlungN(I);
            
            % class branches of right upper lobe
            classedgelobes(RULN, 'RU')
            
            % % Identify mid and lower lobes in right lung
            switch length(RlungN) 
                case 3 % all 3 lobes branch off same node (low res CT)
                    G_copy = G; % make copy and remove upper lobe
                    G_copy = rmnode(G_copy,bfsearch(G_copy,RULN));
                    Gsub = subgraph(G_copy,bfsearch(G_copy,(find(G_copy.Nodes.label == rightN))));
                case 2 % Much more likely
                    % check if following node also belongs to upper lobe
                    upper_ratio = G.Edges.Weight(findedge(G,rightN,RlungN(RlungN ~= RULN)));
                    lower_ratio = G.Edges.Weight(findedge(G,obj.carina_node,rightN));
                    if upper_ratio/lower_ratio < 0.5
                        % non standard upper lobe branching
                        RlungN2 = find(G, RlungN(RlungN ~= RULN));
                        [~, I] = max([obj.Gnode(RlungN2).comz]);
                        RULN2 = RlungN(I);
                        classedgelobes(RULN2, 'RU')
                        Gsub = subgraph(G,bfsearch(G,RlungN2(RlungN2~=RULN2)));
                    else
                    % create subgraph of following nodes
                    Gsub = subgraph(G,bfsearch(G,RlungN(RlungN ~= RULN)));
                    end
                otherwise
                    error('Lobe classification failed to distinguish right lobes, Airways may be incomplete.')
            end

            
            % get list of endpoint nodes
            endpoint = Gsub.Nodes(Gsub.Nodes.ep == 1, :);
            % compute z - y.
            z_minus_y = [endpoint.comz] - [endpoint.comy];
            % identify end node of right middle lobe.
            [~, I] = max(z_minus_y);
            RML_end = endpoint.label(I);
            % identify end node of right lower lobe.
            [~, I] = min(z_minus_y);
            RLL_end = endpoint.label(I);
            % identify bifurcation point of the two end points.
            RML_endpath = flip(shortestpath(G,obj.carina_node, RML_end));
            RLL_endpath = flip(shortestpath(G,obj.carina_node, RLL_end));
            [~, I] = intersect(RML_endpath, RLL_endpath);
            RML_RLLN = RML_endpath(min(I));
            % assign bronchi edges between right node and bifurcating node.
            classedgenonlobes(rightN, RML_RLLN, 'B')
            % identify RML node
            RMLN = RML_endpath(min(I)-1);
            % assign labels to RM edges, not RL
            classedgelobes(RMLN, 'RM')
            
%             % % check for remaining non-lobe branches if they exist.
%             P = shortestpath(G, RightN, RML_RLLN);
%             switch length(P)
%                 case 1 % No further branching from right lung node
%                     % do nothing
%                 case 2 % Further major bronchi to lower lobes
%                     classedgenonlobes(RlungN, RML_RLLN, 'R')
%                 case 3 % potential non-standard branching
%                     upper_ratio = G.Edges.Weight(findedge(G,RightN,P(2)));
%                     lower_ratio = G.Edges.Weight(findedge(G,obj.carina_node,RightN));
%                     if upper_ratio/lower_ratio < 0.5
%                         
%                     end
%             end
            
            % assign remaining labels to RLL
            lobes(cellfun(@isempty,lobes)) = {'RL'};
            
            % add lobe field to glink
            lobe = num2cell(lobes);
            [obj.Glink(:).lobe] = lobe{:};
            
            % save object class
            save(obj)
            
            function classedgelobes(node, label)
                % identify all edges from node outwards
                [~, E] = bfsearch(G, node, 'edgetonew');
                % add predecessing edge
                E = [E; inedges(G, node)];
                lobes(G.Edges.Label(E)) = {label};
            end
            
            function classedgenonlobes(s, t, label)
            [~, ~, E] = shortestpath(G, s, t);
            lobes(G.Edges.Label(E)) = {label};
            end
        end
        
        %% GRAPH NETWORK ANOMOLY ANALYSIS
        
        function [debugseg, debugskel, erroredge] = DebugGraph(obj)
            % TODO: consider moving graph plot to the default plot.
            % First check for multiple 'inedges' to all nodes.
            multiinedgenodes = cell(height(obj.Gdigraph.Edges),1);
            j = 1;
            for i = 1:height(obj.Gdigraph.Nodes)
                eid = inedges(obj.Gdigraph,i);
                if length(eid) > 1
                    multiinedgenodes{j} = eid;
                    j = j + 1;
                end
            end
            multiinedgenodes(j:end)=[];
            % convert cell of different shapes to vector, by Wolfie on
            % stackoverflow.
            
            % 1. Get maximum size of T elements
            %    Pad all elements of T up to maxn values with NaN
            maxn = max(cellfun( @numel, multiinedgenodes ));
            Tpadded = cellfun( @(x) [x; NaN(maxn-numel(x))], multiinedgenodes, 'uni', 0);
            % 2. Convert to array.
            Tpadded = cat(2, Tpadded{:} );
            % 3. Reshape to be one row and remove NaNs
            Trow = reshape( Tpadded.', 1, [] );
            erroredge = Trow(~isnan(Trow));
            
            % return empty debug if no errors found.
            if isempty(erroredge)
                debugseg = [];
                debugskel = [];
                erroredge = [];
                disp('No errors found in skeleton/segmentation')
            else
                warning([int2str(length(multiinedgenodes)) ' errors found in seg/skel. Check output of AirQuant.DebugGraph'])
                % identify edge indices in Glink
                Glink_ind = obj.Gdigraph.Edges.Label(erroredge);
                % get branch labelled segmentation
                branch_seg = ClassifySegmentation(obj);
                % copy over segmentation
                debugseg = double(obj.seg);
                % identify voxels of error branch in seg
                errorvox = ismember(branch_seg, Glink_ind);
                debugseg(errorvox == 1) = 2;
                % identify voxels in skel of 'error' branches
                errorvox = [obj.Glink(Glink_ind).point];
                %copy skel and relabel 'error' branches
                debugskel = double(obj.skel);
                debugskel(errorvox) = 2;
                % also show debug graph plot
%                 figure
%                 G = obj.Gdigraph;
%                 h = plot(G,'EdgeLabel',G.Edges.Label, 'Layout','layered');
%                 h.NodeColor = 'k';
%                 h.EdgeColor = 'k';
%                 highlight(h,'Edges',erroredge,'EdgeColor','r')
                
            end
        end
        
        function [data1, data2] = Remove2nodeinedge(obj, erroredge)
            % identify erroredges that appear more than once 
            % (with same node)
            % TODO: edit to hand >2.
            
            % only assumes that no more than 2 edges can share the same
            % nodes in error.
            erroredge = erroredge';
            
            % get node pairs for each edge
            node_pairs = table2array(obj.Gdigraph.Edges(erroredge,'EndNodes'));
            % identify which ones repeat
            [~,uniqueInd] = unique(node_pairs,'rows','stable');
            duplicateInd = setdiff(1:size(node_pairs,1),uniqueInd);
            % copy one half of the duplicate indices
            duplicateNode = node_pairs(duplicateInd,:);
            duplicateedges = erroredge(duplicateInd);
            
            % Find the other paired duplicate indices.
            paired_edges = cell(length(duplicateedges),1);
            j = 1;
            for i = 1:length(duplicateNode)
                paired_edges{j} = zeros(2,1);
                paired_edges{j}(1,1) = duplicateInd(i);
                [~,~,paired_edges{j}(2,1)] = intersect(duplicateNode(i,:),node_pairs,'rows');
                j = j + 1;
            end
            
            for i = 1:length(paired_edges)
            
            % get skel points of edges and convert to sub
            Glink_ind = obj.Gdigraph.Edges.Label(paired_edges{i});
            edge1 = obj.Glink(Glink_ind(1)).point;
            edge2 = obj.Glink(Glink_ind(2)).point;
            [Y1, X1, Z1] = ind2sub(size(obj.skel), edge1);
            [Y2, X2, Z2] = ind2sub(size(obj.skel), edge2);

            % compute best fit cubic curve
            data1 = [X1', Y1', Z1'];
            data2 = [X2', Y2', Z2'];

            %Generating the spline
            
            end
            % rewrite skeleton with new spline
            
            % overwrite skeleton
            
            % recompute Tree digraph
            
        end
        
        %% UTILITIES
        
        function save(obj)
            save(obj.savename, 'obj', '-v7.3')
        end
        
        function branch_seg = ClassifySegmentation(obj)
            % TODO: consider making this function more robust!
            % label every branch in segmentation to AS branch index.
            
            % Find linear indicies of skeleton
            skel_ind = find(obj.skel == 1);
            classed_skel = zeros(size(skel_ind));
            % for each skeleton branch
            for j = 1:length(obj.Glink)
                % get list of points in branch
                for m = 1:length(obj.Glink(j).point)
                    % put edge number by skeleton index
                    ind = find(skel_ind == obj.Glink(j).point(m));
                    classed_skel(ind) = j;
                end
            end
            
            % get list of everypoint on segmentation
            [XPQ, YPQ, ZPQ] = ind2sub(size(obj.seg),find(obj.seg == 1));
            PQ = [XPQ,YPQ,ZPQ];
            % get list of everypoint on skeleton
            [XP, YP, ZP] = ind2sub(size(obj.seg), skel_ind);
            P = [XP, YP, ZP];
            % find nearest seg point to it on skeleton
            T = delaunayn(P);
            k = dsearchn(P,T,PQ);
            % find that skeleton point's edge assignment
            branch_seg = zeros(size(obj.seg));
            for i = 1:length(PQ)
                branch_seg(PQ(i,1),PQ(i,2),PQ(i,3)) = classed_skel(k(i));
            end
        end
        
        function report = debuggingreport(obj)
            % produce a table showing success of processing for each airway
            % branch.
            report = struct('airway',num2cell(1:length(obj.Glink)));
           
            % check for each airway arclength and FWHM failures
            for i = 1:length(obj.Glink)
                report(i).arclength = ~any(isnan(obj.arclength{i}));
                try
                report(i).FWHM_inner = all(~cellfun(@isempty,obj.FWHMesl{i,1}));
                report(i).FWHM_peak = all(~cellfun(@isempty,obj.FWHMesl{i,2}));
                report(i).FWHM_outer = all(~cellfun(@isempty,obj.FWHMesl{i,3}));
                catch
                end
            end
                        % remove trachea branch. 
            report(obj.trachea_path) = [];
            report = struct2table(report);
        end
                    
        %% HIGH LEVEL METHODS
        function obj = AirwayImageAll(obj)
            % Traverse all airway segments except the trachea.
            disp('Start traversing all airway segments')
            total_branches = length(obj.Glink);
            % check to see if any branches already processed.
            incomplete = cellfun(@isempty, obj.TraversedImage);
            
            for i = 1:length(obj.Glink)
                % skip the trachea or already processed branches
                if i == obj.trachea_path || incomplete(i) == 0
                    continue
                end
                obj = CreateAirwayImage(obj, i);
                disp(['Traversing: Completed ', num2str(i), ' of ', num2str(total_branches)])
            end
            disp('Traversing: Done')
        end
        
        
        function obj = FindFWHMall(obj)
            % Clear all FWHMesl cells
            obj.FWHMesl = cell(length(obj.Glink),3);
            % analyse all airway segments except the trachea.
            disp('Start computing FWHM boundaries of all airway segments')
            total_branches = length(obj.Glink);
            
            % check which branches already traversed
            incomplete = cellfun(@isempty, obj.FWHMesl(:,1));

            for i = 1:length(obj.Glink)
                % skip the trachea
                if i == obj.trachea_path || incomplete(i) == 0
                    continue
                end
                obj = FindAirwayBoundariesFWHM(obj, i);
                disp(['FWHMesl: Completed ', num2str(i), ' of ', num2str(total_branches)])
            end
            save(obj)
            disp('FWHMesl: Done')
        end
        
        %% TRAVERSING AIRWAYS METHODS %%%
        function obj = CreateAirwayImage(obj, link_index)
            % Constructs perpendicular images as if travelling along an
            % airway segment in CT image and Segmentation.
            
            % * Compute whole Spline
            spline = ComputeSpline(obj, link_index);
            spline_para_limit = spline.breaks(end);
            spline_points = 0:obj.spline_sampling_interval:spline_para_limit;
            
            % loop along spline
            % TODO: auto calc 133
            TransAirwayImage = zeros(133,133,length(spline_points));
            TransSegImage = zeros(133,133,length(spline_points));
            arc_length = zeros(length(spline_points),1);
            for i = 1:length(spline_points)
                try
                    % * Compute Normal Vector per spline point
                    [normal, CT_point] = AirQuant.ComputeNormal(spline, spline_points(i));
                    % * Interpolate Perpendicular Slice per spline point
                    [TransAirwayImage(:,:,i), TransSegImage(:,:,i)] = ...
                        InterpolateCT(obj, normal, CT_point);
                    % * Compute real arc_length at this spline point
                    arc_length(i) = Arc_length_to_point(spline_points(i),spline);
                catch
                    % TODO identify why arc_length cannot be calculated in
                    % some cases.
                    warning('Failed to interpolate slice/identify arclength')
                    TransAirwayImage(:,:,i) = zeros(133);
                    TransSegImage(:,:,i) = zeros(133);
                    arc_length(i) = NaN;
                end
            end
            % * Replace NaN entries in images with zero.
            TransAirwayImage(isnan(TransAirwayImage)) = 0;
            TransSegImage(isnan(TransAirwayImage)) = 0;
            % * Save traversed image and arclength for each image
            obj.TraversedImage{link_index, 1} = TransAirwayImage;
            obj.TraversedSeg{link_index, 1} = TransSegImage;
            obj.arclength{link_index, 1} = arc_length;
            % save obj to disk after every branch
            save(obj)
        end
        
        
        function [CT_plane, seg_plane] = InterpolateCT(obj, normal, CT_point)
            % Interpolates a CT plane of the image.
            % Based on original function by Kin Quan 2018
            
            % * Construct spatial physical CT grid
            image_sz = size(obj.CT);
            [x_domain , y_domain , z_domain] = ...
                meshgrid(1:image_sz(2),1:image_sz(1),1:image_sz(3));
            x_domain = x_domain*obj.CTinfo.PixelDimensions(1);
            y_domain = y_domain*obj.CTinfo.PixelDimensions(2);
            z_domain = z_domain*obj.CTinfo.PixelDimensions(3);
            
            % * Get plane grid
            basis_vecs = Orthonormal_basis_with_tangent_vector(normal);
            plane_grid = Grids_coords_for_plane(basis_vecs(:,3),...
                basis_vecs(:,2), CT_point, obj.physical_plane_length,...
                obj.physical_sampling_interval);
            
            % * Execute cubic inperpolation on CT
            plane_intensities = interp3(x_domain,y_domain,z_domain,...
                obj.CT,plane_grid.y(:),plane_grid.x(:),...
                plane_grid.z(:),'cubic');
            
            % * Execute cubic inperpolation on CT
            seg_intensities = interp3(x_domain,y_domain,z_domain,...
                double(obj.seg),plane_grid.y(:),plane_grid.x(:),...
                plane_grid.z(:),'cubic');
            
            % Reshape
            % TODO: Look at what these two lines do.
            plane_length = sqrt(length(plane_grid.y(:)));
            CT_plane = reshape(plane_intensities,...
                [plane_length plane_length]);
            seg_plane = reshape(seg_intensities,...
                [plane_length plane_length]);
        end
        
        
        function spline = ComputeSpline(obj, link_index)
            % Computes a smooth spline of a single graph edge.
            % Based on original function by Kin Quan 2018
            
            % The input is the list of ordered index
            % The output is the smooth spline as a matlab sturct
            
            %Convert into 3d Points
            [x_point, y_point, z_point] = ind2sub(size(obj.CT),obj.Glink(link_index).point);
            
            %Smoothing the data
            voxel_size = obj.CTinfo.PixelDimensions;
            smooth_x = smooth(x_point*voxel_size(1));
            smooth_y = smooth(y_point*voxel_size(2));
            smooth_z = smooth(z_point*voxel_size(3));
            
            %Complete smooth data
            smooth_data_points = [smooth_x smooth_y smooth_z]';
            
            %Generating the spline
            spline = cscvn(smooth_data_points);
        end
        
        %% MEASUREMENT METHODS
        function obj = FindAirwayBoundariesFWHM(obj, link_index)
            %Based on function by Kin Quan 2018 that is based on Kiraly06
            
            slices_size = size(obj.TraversedImage{link_index, 1}, 3);
            
            % Prepping the outputs
            raycast_FWHMl = cell(slices_size, 1);
            raycast_FWHMp = cell(slices_size, 1);
            raycast_FWHMr = cell(slices_size, 1);
            
            % For every traversed slice
            for k = 1:slices_size
                try
                    % * Compute airway centre
                    % Check that airway centre is slice centre
                    center = ...
                        Return_centre_pt_image(...
                        obj.TraversedSeg{link_index, 1}(:,:,k));
                    
                    % Recompute new centre if necessary
                    [centre_ind , new_centre] =  ...
                        Check_centre_with_segmentation(obj.TraversedImage{link_index, 1}(:,:,k), ...
                        obj.TraversedSeg{link_index, 1}(:,:,k));
                    if ~centre_ind
                        center = fliplr(new_centre);
                    end
                    
                    % * Raycast
                    [CT_rays, seg_rays, coords]= Raycast(obj, ...
                        obj.TraversedImage{link_index, 1}(:,:,k), ...
                        obj.TraversedSeg{link_index, 1}(:,:,k), center);
                    
                    % * Compute FWHM
                    [FWHMl, FWHMp, FWHMr] = AirQuant.computeFWHM(CT_rays,...
                        seg_rays, coords);
                    
                    % * Compute Ellipses
                    FWHMl_ellipse = ComputeEllipses(obj, FWHMl);
                    FWHMp_ellipse = ComputeEllipses(obj, FWHMp);
                    FWHMr_ellipse = ComputeEllipses(obj, FWHMr);
                    
                    % * Record, catch incase a slice fails.
                    raycast_FWHMl{k,1} = FWHMl_ellipse;
                    raycast_FWHMp{k,1} = FWHMp_ellipse;
                    raycast_FWHMr{k,1} = FWHMr_ellipse;
                catch
                    % segmentation exceeds interpolated slice therefore no
                    % measurement recorded.
                    raycast_FWHMl{k,1} = NaN;
                    raycast_FWHMp{k,1} = NaN;
                    raycast_FWHMr{k,1} = NaN;
                end
            end
            obj.FWHMesl{link_index, 1} = raycast_FWHMl;
            obj.FWHMesl{link_index, 2} = raycast_FWHMp;
            obj.FWHMesl{link_index, 3} = raycast_FWHMr;
        end
        
        function [CT_rays, seg_rays, coords] = Raycast(obj, interpslice, interpseg, center)
            % * Compute Rays
            % Getting the range limit of the ray which will be the shortest
            % distance from the centre to the bounadry Need to find the
            % limits of the raw
            image_size = size(interpslice);
            limit_row = abs(image_size(1) - center(1));
            limit_col = abs(image_size(2) - center(2));
            ray_length_limit = min(limit_row, limit_col);
            
            %Compute rays in polar coords
            ray_angle_interval = 2*pi/obj.num_rays;
            radial = 0:obj.ray_interval:ray_length_limit;
            theata = 0:ray_angle_interval:2*pi;
            
            %Convert rays to cartesian coords
            x_component = radial'*cos(theata) + center(1);
            y_component = radial'*sin(theata) + center(2);
            
            % * Cast rays
            interpslice = double(interpslice);
            
            CT_rays = interp2(interpslice, x_component(:),...
                y_component(:));
            
            seg_rays = interp2(interpseg, x_component(:),...
                y_component(:));
            
            %Need to reshape
            CT_rays = ...
                reshape(CT_rays,[size(y_component,1) size(y_component,2)]);
            seg_rays = ...
                reshape(seg_rays,[size(y_component,1) size(y_component,2)]);
            coords = cat(3, x_component, y_component);
        end
        
        function elliptical_sturct = ComputeEllipses(obj, FWHM_points)
            %Perform the Ellipitcal fitting
            %The output will be sturct containing the major and minor lengths as well
            %as all the stop points
            elliptical_sturct = struct;
            elliptical_sturct.x_points = FWHM_points(:,1);
            elliptical_sturct.y_points = FWHM_points(:,2);
            % Compute fitting
            elliptical_sturct.elliptical_info = ...
                Elliptical_fitting(elliptical_sturct.x_points, ...
                elliptical_sturct.y_points);
            elliptical_sturct.area = ...
                elliptical_sturct.elliptical_info(3)*...
                elliptical_sturct.elliptical_info(4)*...
                pi*obj.physical_sampling_interval.^2;
        end
     
        %% ANALYSIS METHODS
%         function obj = ComputeTaperValues(obj)
%             for i = 1:length(obj.TraversedImage)
%                 cum_area = [];
%                 for j = 1:length(obj.arclength{i,1})
%                     try
%                         cum_area = [cum_area; obj.FWHMesl{i, 1}{j, 1}.area];
%                     catch
%                         cum_area = [cum_area; NaN];
%                     end
%                 end
%                 obj.Specs(i).FWHMl_logtaper = ...
%                     AirQuant.ComputeTaperGrad(...
%                     obj.arclength{i, 1}, cum_area);
%             end
%         end
        
        function [logtapergrad, cum_arclength, cum_area, edgepath] = ConstructTaperPath(obj, terminal_node_idx, type)
            if nargin < 3
                type = 'inner';
            end
            switch type
                case 'inner'
                    typeidx = 1;
                case 'peak'
                    typeidx = 2;
                case 'outer'
                    typeidx = 3;
                otherwise
                    error(['Unknown type: ', type, 'please use ''inner'', ''peak'' or ''outer''.'])
            end
            % TODO: remove trachea node?
            G = digraph(obj);
            [path,~,~] = shortestpath(G, obj.carina_node, terminal_node_idx);
            
            % construct edgepath from graph2skel NB: MATLAB graph object
            % does not work.
            nodelist = [[obj.Glink(:).n1]', [obj.Glink(:).n2]'];
            edgepath = zeros(length(path)-1, 1);
            for i = 1:length(path)-1
                nodepair = [path(i), path(i+1)];
                [~,edgepath(i)] = ismember(nodepair, nodelist,'rows');
            end
            cum_arclength = 0;
            cum_area = [];
            k = 1;
            for q = 1:length(edgepath)
                i = edgepath(q);
                % skip the first of each link except the first one
                if k == 1
                    try
                        cum_area = [cum_area; obj.FWHMesl{i, typeidx}{1, 1}.area];
                        k=0;
                    catch
                        cum_area = [cum_area; NaN];
                        k=0;
                    end
                end
                % Appending arclengths
                max_current_arclength = max(cum_arclength);
                current_arclength_array = max_current_arclength+obj.arclength{i,1}(2:end);
                cum_arclength = [cum_arclength; current_arclength_array];
                
                for j = 2:length(obj.arclength{i,1})
                    try
                        cum_area = [cum_area; obj.FWHMesl{i, typeidx}{j, 1}.area];
                    catch
                        cum_area = [cum_area; NaN];
                    end
                end
            end
            logtapergrad = AirQuant.ComputeTaperGrad(cum_arclength, cum_area);
            
            % TODO: consider using graph and edge highlight for another
            % applications.
            %h = plot(G);
            %highlight(h,'Edges',edgepath,'EdgeColor','r','LineWidth',1.5)
        end
        
        
        function terminalnodelist = ListTerminalNodes(obj)
            terminalnodelist = find([obj.Gnode(:).ep] == 1);
            % remove trachea node
            terminalnodelist(terminalnodelist ==  obj.trachea_node) = [];
        end
        
        
        function AllTaperResults = ComputeTaperAll(obj)
            % get list of terminal branches
            terminallinklist = ListTerminalNodes(obj);
            % construct structure to save analysis results
            AllTaperResults = struct('terminalnode',num2cell(terminallinklist));
            % compute taper results for each terminal branch
            for i = 1:length(AllTaperResults)
                node_idx = AllTaperResults(i).terminalnode;
                
                [logtapergrad_inner, cum_arclength, cum_area_inner, edgepath] = ...
                    ConstructTaperPath(obj, node_idx, 'inner');
                [logtapergrad_peak, ~, cum_area_peak, ~] = ...
                    ConstructTaperPath(obj, node_idx, 'peak');
                [logtapergrad_outer, ~, cum_area_outer, ~] = ...
                    ConstructTaperPath(obj, node_idx, 'outer');
                AllTaperResults(i).edgepath = edgepath;
                AllTaperResults(i).arclength = cum_arclength;
                AllTaperResults(i).area_inner = cum_area_inner;
                AllTaperResults(i).area_peak = cum_area_peak;
                AllTaperResults(i).area_outer = cum_area_outer;
                AllTaperResults(i).logtapergrad_inner = logtapergrad_inner;
                AllTaperResults(i).logtapergrad_peak = logtapergrad_peak;
                AllTaperResults(i).logtapergrad_outer = logtapergrad_outer;
                
                try % only add lobe information if it exists
                    AllTaperResults(i).lobe = obj.Glink(obj.Gnode(node_idx).links).lobe;
                catch
                end
            end
            if isfield(AllTaperResults, 'lobe')
                AllTaperResults = sortrows(struct2table(AllTaperResults), 'lobe');
            else
                AllTaperResults = struct2table(AllTaperResults);
            end
        end
        
        
        function PlotTaperResults(obj, terminal_node_idx, type)
            if nargin < 3
                type = 'other';
            end
            switch type
                case 'inner'
                    logtapergrad = plottaperunderfunc(obj, terminal_node_idx, type);
                    typetxt = 'Inner lumen';
                case 'peak'
                    logtapergrad = plottaperunderfunc(obj, terminal_node_idx, type);
                    typetxt = 'Peak wall';
                case 'outer'
                    logtapergrad = plottaperunderfunc(obj, terminal_node_idx, type);
                    typetxt = 'Outer wall';
                otherwise
                    subplot(3,1,1)
                    PlotTaperResults(obj, terminal_node_idx, 'inner')
                    subplot(3,1,2)
                    PlotTaperResults(obj, terminal_node_idx, 'peak')      
                    subplot(3,1,3)
                    PlotTaperResults(obj, terminal_node_idx, 'outer')
                    return % end function
            end
            ylabel('Area (mm^2)')
            txt = sprintf('%s log taper graph; Terminal-node: %d; LogTaperGrad: %0.3g',typetxt, terminal_node_idx,logtapergrad);
            title(txt)
            
            function logtapergrad = plottaperunderfunc(obj, terminal_node_idx, type)
            [logtapergrad, cum_arclength, cum_area, ~] = ConstructTaperPath(obj, terminal_node_idx, type);
            plot(cum_arclength, cum_area, 'k.')
            hold on
            [~, displacement] = AirQuant.ComputeTaperGrad(cum_arclength, cum_area);
            logcurve = exp(-logtapergrad*cum_arclength + displacement);
            plot(cum_arclength, logcurve,'-r')
            xlabel('Arc-length (mm)')
            legend('Measured', 'Log curve fit')
        end
        end
        

        
        %% Graph Methods
        function G = digraph(obj)
            % TODO: REDUNDANT FUNCTION?
            % compute edge weights
            weights = zeros(length(obj.Glink), 1);
            for i = 1:length(obj.Glink)
                weights(i) = length(obj.Glink(i).point);
            end
            % add edges
            G = digraph([obj.Glink(:).n1],[obj.Glink(:).n2], weights);
        end
        
        %% VISUALISATION METHODS
        function h = plot(obj, type)
            % Default plot is a graph network representation. Optional input is to
            % provide a list of edge labels indexed by the Glink property.
            
            % Only show graph for airways from carina to distal.
            G = obj.Gdigraph;
            
            %     trachea_edges = find(G.Edges.Label == obj.trachea_path);
            %     G = rmedge(G, trachea_edges);
            %     G = rmnode(G, find(indegree(G)==0 & outdegree(G)==0));
            
            if nargin > 1
                labeltype = type;
            else
                % default edge label is Glink index.
                labeltype = 'index';
            end
            
            switch labeltype
                case 'index' % default
                    edgelabels = G.Edges.Label;
                case 'lobes'
                    try
                        lobes = [obj.Glink(:).lobe];
                    catch
                        error('Need to run ComputeAirwayLobes first')
                    end
                case {'generation','gen'}
                    gens = [obj.Glink(:).generation];
                    edgelabels = gens(G.Edges.Label);
                otherwise
                    warning('Unexpected plot type. Resorting to default type.')
                    edgelabels = G.Edges.Label;
            end
            
            h = plot(G,'EdgeLabel',edgelabels, 'Layout', 'layered');
            h.NodeColor = 'r';
            h.EdgeColor = 'k';
            % check for error edges
            [~,~,erroredge] = DebugGraph(obj);
            if ~isempty(erroredge)
                highlight(h,'Edges',erroredge,'EdgeColor','r')
            end
        end
        
        function h = plotgenlabels(obj)
            gens = [AS.Glink(:).generation];
            genslabels = gens(obj.Gdigraph.Edges.Label);
            h = plot(obj, genslabels);
        end
        
        
        function PlotTree(obj)
            % Plot the airway tree with nodes and links
            % Original Function by Ashkan Pakzad on 27th July 2019.
            
            isosurface(obj.skel);
            alpha(0.7)
            hold on
            
            % edges
            ind = zeros(length(obj.Glink), 1);
            for i = 1:length(obj.Glink)
                ind(i) = obj.Glink(i).point(ceil(end/2));
            end
            [Y, X, Z] = ind2sub(size(obj.seg),ind);
            nums_link = string(1:length(obj.Glink));
            %plot3(X,Y,Z, 'b.', 'MarkerFaceColor', 'none');
            text(X+1,Y+1,Z+1, nums_link, 'Color', [0, 0.3, 0])
            
            % nodes
            X_node = [obj.Gnode.comy];
            Y_node = [obj.Gnode.comx];
            Z_node = [obj.Gnode.comz];
            nums_node = string(1:length(obj.Gnode));
            plot3(X_node,Y_node,Z_node, 'r.', 'MarkerSize', 18, 'Color', 'r');
            text(X_node+1,Y_node+1,Z_node+1, nums_node, 'Color', [0.8, 0, 0])
            
            %axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            view(80,0)
            axis vis3d
            
        end
        
        function PlotMap3D(obj, mode)
            % mode = 'TaperGradient', 'generation', 'lobes'
            axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            
            % generating the color data
            cdata = zeros(size(obj.seg));
            branch_seg = ClassifySegmentation(obj);
            switch mode 
                case 'TaperGradient'
                    % TODO: rewrite this bit....
                    for i = 1:length(obj.Specs)
                        cdata(branch_seg == i) = obj.Specs(i).FWHMl_logtaper*-1;
                    end
                case 'Generation'
                    for i = 1:length(obj.Glink)
                        cdata(branch_seg == i) = obj.Glink(i).generation;
                    end
                    clims = [0 max(cdata(:))];
%                     colorbarstring = 'Generation Number';
                case 'Lobe'
                    % convert lobe id to number
                    lobeid = {'B','RU','RM','RL','LU','lin','LL'};
                    for i = 1:length(obj.Glink)
                        cdata(branch_seg == i) = find(strcmp(lobeid, obj.Glink(i).lobe))-1;
                    end
                    clims = [0 max(cdata(:))+1];
%                     colorbarstring = 'Lobe';
%                     colourshow = clims(1):clims(2);
%                     colourlabels = lobeid;
                    
                otherwise
                    error('Choose appropiate mode.')
            end
            % producing segmentation 3d object
            p = patch(isosurface(obj.seg));
            % map colour indices to 3d object vertices
            isocolors(cdata, p);
            % % Reassign colors to meaningful values
            Vertcolour = p.FaceVertexCData;
            vcolours = uniquetol(Vertcolour);
            colourind = unique(cdata);
            for i = 1:length(vcolours)
                Vertcolour(Vertcolour==vcolours(i)) = colourind(i);
            end
            % overwrite vertices colour with meaningful values
            p.FaceVertexCData = Vertcolour;
            
            p.FaceColor = 'interp';
            p.EdgeColor = 'none';
            map = linspecer(max(cdata(:))+2);
            colormap(map)
%             c = colorbar('Ticks', colourshow, 'TickLabels', colourlabels);
%             c.Label.String = colorbarstring;
            caxis(clims)
            view(80,0)
            axis vis3d
        end
        
        function PlotAirway3(obj, link_index)
            % Plot resampled airway slices overlayed with FWHMesl ray cast
            % points and fitted ellipse
            f = figure('Position',  [100, 100, 850, 600]);
            slide = 1;
            PlotAirway(obj, link_index, slide)
            numSteps = size(obj.TraversedImage{link_index,1}, 3);
            
            b = uicontrol('Parent',f,'Style','slider','Position',[50,10,750,23],...
                'value',slide, 'min',1, 'max',numSteps, 'SliderStep', [1/(numSteps-1) , 1/(numSteps-1)]);
            bgcolor = f.Color;
            uicontrol('Parent',f,'Style','text','Position',[25,10,23,23],...
                'String', '1','BackgroundColor',bgcolor);
            uicontrol('Parent',f,'Style','text','Position',[800,10,23,23],...
                'String',numSteps,'BackgroundColor',bgcolor);
            
            b.Callback = @sliderselect;
            
            function sliderselect(src,event)
                val=round(b.Value);
                PlotAirway(obj, link_index, val);
            end
            
        end
        
        function PlotAirway(obj, link_index, slide)
            % display image
            imagesc(obj.TraversedImage{link_index, 1}(:,:,slide))
            hold on
            colormap gray
            
            try % try block incase FWHMesl has not been executed.
                % plot ray cast results
                
                plot(obj.FWHMesl{link_index, 1}{slide, 1}.x_points, obj.FWHMesl{link_index, 1}{slide, 1}.y_points,'r.')
                plot(obj.FWHMesl{link_index, 2}{slide, 1}.x_points, obj.FWHMesl{link_index, 2}{slide, 1}.y_points,'c.')
                plot(obj.FWHMesl{link_index, 3}{slide, 1}.x_points, obj.FWHMesl{link_index, 3}{slide, 1}.y_points,'y.')
                
                % plot ellipse fitting
                ellipse(obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(4),...
                    obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(1),...
                    obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(2),'m');
                
                ellipse(obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(4),...
                    obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(1),...
                    obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(2),'b');
                
                ellipse(obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(4),...
                    obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(1),...
                    obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(2),'y');
                %TODO: set third colour more appropiately
                
            catch
                % warning('No FWHMesl data, showing slices without elliptical information.')
            end
            
            % display area measurements
            % TODO: is this needed?
            %dim = [.15 .85 .24 .05];
            %a = annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','y');
            a = rectangle('Position',[0,0,133,10],'FaceColor','y','LineWidth',2);
            ax = gca;
            try
                text(ax, 1,5,sprintf('Arc Length = %4.1f mm; Inner area = %4.2f mm^2; Peak area = %4.2f mm^2; Outer area = %4.2f mm^2; %3.0i of %3.0i', ...
                    obj.arclength{link_index, 1}(slide), obj.FWHMesl{link_index, 1}{slide, 1}.area, obj.FWHMesl{link_index, 2}{slide, 1}.area ,...
                    obj.FWHMesl{link_index, 3}{slide, 1}.area, slide, size(obj.TraversedImage{link_index, 1},3)));
            catch
                text(ax, 1,5,sprintf('Arc Length = %4.1f mm; %3.0i of %3.0i', ...
                    obj.arclength{link_index, 1}(slide), slide, size(obj.TraversedImage{link_index, 1},3)));
            end
        end
        
    end
    %% STATIC METHODS
    methods (Static)
        function [normal, CT_point] = ComputeNormal(spline, point)
            % Based on original function by Kin Quan 2018
            % * interperate real point on spline
            CT_point = fnval(spline, point);
            % * get tangent of point along spline
            % differentiate along spline to get gradient
            spline_1diff = fnder(spline,1);
            tangent_vec = fnval(spline_1diff,point);
            normal = tangent_vec/norm(tangent_vec,2);
        end
        
        
        function [logtapergradient, displacement] = ComputeTaperGrad(arclength, area)
            % identify NaN data points
            idx = isnan(area);
            % compute logtapergrad, ignoring nan values
            p_coeff = polyfit(arclength(~idx),log(area(~idx)),1);
            % positive if thinning, therefore multiply by -1.
            logtapergradient = p_coeff(1) * -1; % logtapergradient
            if nargout > 1
                displacement = p_coeff(2); % displacement
            end
        end
        
        
        function [FWHMl, FWHMp, FWHMr] = computeFWHM(CT_rays, seg_rays, coords)
            %This is to perfrom the ray casting measurents - the input is in a sturct
            % as
            
            % The basis of the code is base on the description - Virtual Bronchoscopy for Quantitative Airway Analysis
            % by A P Kiraly et al.
            
            %Need to find the length of the ray - this need a loop
            number_of_rays = size(CT_rays,2);
            FWHMl_x = [];
            FWHMl_y = [];
            FWHMp_x = [];
            FWHMp_y = [];
            FWHMr_x = [];
            FWHMr_y = [];
            
            %performing a loop
            for ray = 1:number_of_rays
                
                CT_profile = CT_rays(:,ray);
                seg_profile = seg_rays(:,ray);
                
                % * Find the edge of the interpolated segmentation.
                if seg_profile(1) < 0.5
                    continue % skip if seg edge does not exist
%                 elseif seg_profile(length(seg_profile)) > 0.5
%                     % set to slice edge if segmentation exceeds slice
%                     seg_half = length(ind_ray);
                else
                    % Identify the last vaule that is above the 0.5
                    ind_ray = (seg_profile < 0.5);
                    seg_half = find(ind_ray,1);
                end
                
                % * Find FWHM peak
                [max_int_array , max_location_array] = ...
                    findpeaks(CT_profile);
                
                if isempty(max_int_array)
                    continue % skip ray if peak does not exist
                end
                
                % identify the peak closest to seg half point
                index_diff = abs(max_location_array - seg_half);
                [~,closest_max] = min(index_diff);
                % ensure the point is unique
                FWHMp = max_location_array(closest_max(1));
                
                % * Compute FWHM on left side
                % loop through profile from FWHM peak to centre
                % find first minima from right to left.
                for i = FWHMp:-1:2
                    if ~(CT_profile(i) >= CT_profile(i - 1))
                        break
                    end
                    i = 1; % incase inner is at beginning of profile
                end
                FWHMi = i; % inner FWHM curve
                
                threshold_int_left = ...
                    (CT_profile(FWHMp) + CT_profile(FWHMi))/2;
                FWHMl = Finding_midpoint_stop(CT_profile,threshold_int_left,FWHMi,FWHMp);
                
                % * Compute FWHM on right side
                % loop through profile from FWHM peak to distal
                % find first minima from left to right.
                for i = FWHMp:length(CT_profile)-1
                    if CT_profile(i) <= CT_profile(i + 1)
                        break
                    end
                    i = length(CT_profile);
                end
                FWHMo = i;  % outer FWHM curve
                % TODO: move threshold calculation into function.
                threshold_int_right = ...
                    (CT_profile(FWHMp) + CT_profile(FWHMo))/2;
                FWHMr = Finding_midpoint_stop_right(CT_profile,threshold_int_right,FWHMo,FWHMp);
                
                % TODO: add same anomaly detection as before.
                
                % concat points together
                FWHMl_x = cat(1,FWHMl_x,coords(FWHMl, ray, 1));
                FWHMl_y = cat(1,FWHMl_y,coords(FWHMl, ray, 2));
                FWHMp_x = cat(1,FWHMp_x,coords(FWHMp, ray, 1));
                FWHMp_y = cat(1,FWHMp_y,coords(FWHMp, ray, 2));
                FWHMr_x = cat(1,FWHMr_x,coords(FWHMr, ray, 1));
                FWHMr_y = cat(1,FWHMr_y,coords(FWHMr, ray, 2));
            end
            
            FWHMl = cat(2, FWHMl_x, FWHMl_y);
            FWHMp = cat(2, FWHMp_x, FWHMp_y);
            FWHMr = cat(2, FWHMr_x, FWHMr_y);
        end
    end
end