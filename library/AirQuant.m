classdef AirQuant < handle % handle class
    properties
        %%% Volumes and Metadata
        CT % CT image
        CTinfo % CT metadata
        seg % binary airway segmentation
        skel % skeleton based on segementation
        savename % filename to load/save
        %%% CT resampling params
        max_plane_sz = 40;% max interpolated slice size
        plane_sample_sz % interpolated slice pixel size
        spline_sample_sz % mm interval to sample along branch arclength
        plane_scaling_sz = 5; % scale airway diameter approx measurement.
        min_tube_sz % smallest measurable lumen Diameter.
        %%% Ray params
        num_rays = 50; % check methods
        ray_interval = 0.2; % check methods

    end
    properties (SetAccess = private)
        %%% Graph Properties
        Gadj % undirected Graph Adjacency matrix
        Gnode % Graph node info
        Glink % Graph edge info
        Gdigraph % digraph object
        trachea_path % edges that form a connect subgraph above the carina
        carina_node % node that corresponds to the carina
        skel_points % list of skeleton points
        %%% Resampled image slices along graph paths/airway segments.
        Dmap % distance transform of seg
        Splines % spline data of branches
        TraversedImage % perpendicular CT slices of all airways
        TraversedSeg % perpendicular segmentation slices of all airways
        arclength % corresponding arclength measurement traversed slices
        FWHMesl % FWHMesl algorithm results for every airway {inner, peak, outer}
        Specs % Store of end metrics
    end
    
    methods
        %% INITIALISATION METHODS
        % Methods used to call and make AQ object before any further
        % processing takes place.
        function obj = AirQuant(CTimage, CTinfo, segimage, skel, savename)
            % Initialise the AirQuant class object.
            % if using default settings, do not provide params argument.
            
            % Can provide just file name, if analyses previously complete.
            if nargin == 1
                savename = CTimage;
                obj.savename = savename;
            end
            
            if isfile(savename)
                disp(['Found case stored at ', savename, ' loading ...'])
                load(savename)
                obj.savename = savename;
            else 
                disp(['New case to be saved at ', savename, ' saving...'])
                obj.CTinfo = CTinfo;
                obj.CT = reorientvolume(CTimage, obj.CTinfo);                
                % TODO: consider preprocess segmentation to keep largest
                % connected component.
                % ensure no holes in segmentation
                obj.seg = reorientvolume(imfill(segimage,'holes'), obj.CTinfo);
                % Set Resampling parameters and limits
                measure_limit = floor((min(obj.CTinfo.PixelDimensions)/2)*10)/10;
                obj.plane_sample_sz = measure_limit;
                obj.spline_sample_sz = measure_limit;
                obj.min_tube_sz = 3*max(obj.CTinfo.PixelDimensions);
                % graph airway skeleton
                obj = GenerateSkel(obj,skel);
                % Convert into digraph
                obj = AirwayDigraph(obj);
                % Identify Carina
                obj = FindCarina(obj);
                % Identify paths that belong to trachea
                obj = FindTracheaPaths(obj);
                % classify airway generations
                obj = ComputeAirwayGen(obj);
                % Compute distance transform
                obj.Dmap = bwdist(~obj.seg);
                % set up empty cell for traversed CT and segmentation
                obj.Splines = cell(length(obj.Glink),2);
                obj.TraversedImage = cell(length(obj.Glink),1);
                obj.TraversedSeg = cell(length(obj.Glink),1);
                % set up empty specs doubles
                obj.arclength = cell(length(obj.Glink),1);
                obj.Specs = struct();
                % set up empty cell for recording raycast/fwhmesl method
                obj.FWHMesl = cell(length(obj.Glink),3);
                try
                    ComputeAirwayLobes(obj)
                catch
                    warning('Airway classification by Lobe unsuccessful. Lobe dependent functions will not work.')
                end
                % save class
                obj.savename = savename;
                save(obj)
                disp('New case AirQuant Object successfully initialised.')
            end
        end
        
        function obj = GenerateSkel(obj, skel)
            % Generate the airway skeleton
            % if skeleton not provided, generate one.
            if isempty(skel)
                obj.skel = Skeleton3D(obj.seg);
            else
                obj.skel = reorientvolume(skel, obj.CTinfo);
            end
            % create graph from skeleton.
            [obj.Gadj,obj.Gnode,obj.Glink] = Skel2Graph3D(obj.skel,0);
            
        end
       
        
        function obj = AirwayDigraph(obj)
            % Converts the output from skel2graph into a digraph tree
            % network, such that there are no loops and edges point from
            % the trachea distally outwards.
            
            % Identify the Trachea node. FindFWHMall
            % Assumes trachea fully segmented and towards greater Z.
            % Smoothen
            [~, trachea_node] = max([obj.Gnode.comz]);
            
            % Create digraph with edges in both directions, loop through
            % and remove if found later in the BF search.
            G = digraph(obj.Gadj);
            % BF search from carina node to distal.
            node_discovery = bfsearch(G,trachea_node);
            %%% reorder nodes by bfs
            G = reordernodes(G, node_discovery);
            % reorder in gnodes and glinks
            obj.Gnode = obj.Gnode(node_discovery);
            for iinode = 1:length(obj.Gnode)
            conns = obj.Gnode(iinode).conn;
                for iiconn = 1:length(conns)
                    obj.Gnode(iinode).conn(iiconn) = find(node_discovery == conns(iiconn));
                end
            end
            for iilink = 1:length(obj.Glink)
                obj.Glink(iilink).n1 = find(node_discovery == obj.Glink(iilink).n1);
                obj.Glink(iilink).n2 = find(node_discovery == obj.Glink(iilink).n2);
            end
            % NB: Trachea node now always  = 1.
                          
            % organise into peripheral facing digraph
            % half of the edges will be removed
            removal = zeros(height(G.Edges)/2 , 1);
            j = 1;
            for i = 1:height(G.Edges)
                % edges should be connected by smaller to larger node.
                if (G.Edges.EndNodes(i,1) - G.Edges.EndNodes(i,2)) > 0
                    % record if not central-->distal
                    removal(j) = i;
                    j = j + 1;
                end
            end
            % remove recorded edges at end in opposite direction.
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
                    continue
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
            for i = 1:length(obj.trachea_path)
                lobes{obj.trachea_path(i)} = 'B';
            end
            
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

            % assign branch from left lobe divider to upper itk node
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
            classedgelobes(LN, 'LUlin')
                        
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
                        RlungN2 = successors(G, RlungN(RlungN ~= RULN));
                        [~, I] = max([obj.Gnode(RlungN2).comz]);
                        RULN2 = RlungN2(I);
                        classedgelobes(RULN2, 'RU')
                        Gsub = subgraph(G,bfsearch(G,RlungN2(RlungN2~=RULN2)));
                    else
                        % create subgraph of following nodes
                        Gsub = subgraph(G,bfsearch(G,RlungN(RlungN ~= RULN)));
                    end
                otherwise
                    error('Lobe classification failed to distinguish right lobes, airway segmentation may be incomplete.')
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
            
            % assign remaining labels to RLL
            lobes(cellfun(@isempty,lobes)) = {'RL'};
            
            % add lobe field to glink
            lobe = num2cell(lobes);
            [obj.Glink(:).lobe] = lobe{:};
           
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
        % Methods that analyse the airway structural graph checking for
        % anomalies such as loops and reporting them.
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
        % Useful functions.
        function save(obj)
            save(obj.savename, 'obj', '-v7.3')
        end
        
        function obj = ComputeSkelPoints(obj)
            [XP, YP, ZP] = ind2sub(size(obj.seg), find(obj.skel == 1));
            obj.skel_points = [XP, YP, ZP]; % list of skel points
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
            if isempty(obj.skel_points)
                ComputeSkelPoints(obj)
            end
            P = obj.skel_points;
            % find nearest seg point to it on skeleton
            T = delaunayn(P);
            k = dsearchn(P,T,PQ);
            % find that skeleton point's edge assignment
            branch_seg = zeros(size(obj.seg));
            for i = 1:length(PQ)
                branch_seg(PQ(i,1),PQ(i,2),PQ(i,3)) = classed_skel(k(i));
            end
        end
        
        function out = SuccessReport(obj)
            % produce a table showing success of processing for each airway
            % branch.
            report = struct('airway',num2cell(1:length(obj.Glink)));
            notanyall = @(x) ~any(x,'all');
            % check for each airway arclength and FWHM failures
            for i = 1:length(obj.Glink)
                if i == obj.trachea_path
                    continue
                end
                % no nan or empty entries arclength
                report(i).arclength = ~any(isnan(obj.arclength{i})) & ~any(isempty(obj.arclength{i}));
                % no nan or empty entries interpolated ct
                intepolatedCT_nonan = all(cellfun(notanyall, cellfun(@isnan,obj.TraversedImage{i,1},'UniformOutput',false)));
                interpolatedCT_noempty = all(cellfun(notanyall, cellfun(@isempty,obj.TraversedImage{i,1},'UniformOutput',false)));
                report(i).InterpolatedCT = intepolatedCT_nonan & interpolatedCT_noempty;
                try
                % nans acceptable, none empty
                report(i).FWHM_inner = ~isempty(obj.FWHMesl{i,1});
                report(i).FWHM_peak = ~isempty(obj.FWHMesl{i,2});
                report(i).FWHM_outer = ~isempty(obj.FWHMesl{i,3});
                catch
                end
            end
                        % remove trachea branch. 
            report(obj.trachea_path) = [];
            report = struct2table(report);
            disp('Debugging report. Success of executing each stage of the algorithm for each airway.')
            disp(report)
            if nargout > 0 
                out = report;
            end
        end
                  
        function [out, h] = AirwayCounts(obj, type, showfig)
            % produce a table and barchart showing 
            % number of airways per generation
            % number of airways per generation per lobe
            % type is either generation or lobe
            % showfig is a flag for plotting true by default
            
            if nargin < 3
                showfig = 1;
            end
            
            % get generation index for each branch
            awygen = [obj.Glink(:).generation]';
            % get 1:max generations
            generations = unique(awygen);
            % set up cell for each gen
            rownames = compose('%d', generations);
            rownames = cellstr(string(rownames));
            
            % count for all gens
            [gencount, ~] = histcounts(awygen, [generations; generations(end)+1]);
            
            switch type
                case 'generation'
                    % Generate count for each generation used bin edges
                    % method
                    if showfig == 1
                        % show figure result as bar chart
                        bar(0:length(gencount)-1, gencount)
                        title('Number of airways per generation')
                        xlabel('Generation')
                        ylabel('count')
                        grid on
                    end
                    
                    reporttable = table(gencount', 'RowNames', rownames, 'VariableNames', {'NGenerations'});
                    disp('Number of airways per generation')
                    disp(reporttable)
                    
                case 'lobe'
                    % also get lobe index of each branch
                    awylobes = [obj.Glink(:).lobe]';
                    % get list of lobe labels and make table
                    varNames = AirQuant.LobeLabels();

                    % loop though each lobe label and generate histogram
                    % count for each generation.
                    lobecount = zeros(length(varNames), length(generations));
                    for iilobe = 1:length(varNames)
                        currentgens = awygen(strcmp(awylobes,varNames{1,iilobe}));
                        [lobecount(iilobe, :), ~] = ...
                            histcounts(currentgens, [generations; ...
                            generations(end)+1]);
                    end
                    reporttable = array2table(lobecount');
                    reporttable.Properties.VariableNames = varNames;
                    reporttable.Properties.RowNames = rownames;
                    disp('Number of airways per generation per lobe')
                    disp(reporttable)
                    if showfig == 1
                        % show figure result as bar chart
                        iilobe = 0;
                        for plotind = [1,3,5,2,4,6]
                            iilobe = iilobe + 1;
                            subplot(3,2,plotind);
                            bar(0:length(gencount)-1, lobecount(iilobe,:));
                            title(varNames(iilobe))
                            xlabel('Generation')
                            ylabel('count')
                            grid on
                        end
                        f = gcf;
                        allax = findall(f,'type','axes');
                        linkaxes(allax,'xy');
                    end
                    
                otherwise
                    error('Choose type: generation or lobe')
            end
            
            if nargout > 0
                out = reporttable;
            end
            if nargout > 1
                h = 1; % change to plot later
            end
        end
        
        %% HIGH LEVEL METHODS
        % methods that package lower level methods, often to apply analysis
        % to all airways rather than just individual airways.
        function obj = AirwayImageAll(obj)
            % Traverse all airway segments except the trachea.
            disp('Start traversing all airway segments')
            total_branches = length(obj.Glink);
            % check to see if any branches already processed.
            incomplete = cellfun(@isempty, obj.TraversedImage);
            
            for i = 1:length(obj.Glink)
                % skip the trachea or already processed branches
                if i == obj.trachea_path || incomplete(i) == 0
                    disp(['Traversing: ', num2str(i), ' trachea skipped or already complete'])
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
        
        %% SPLINE METHODS
        % Technical: These methods compute the central airway spline.
        function ComputeSpline(obj, link_index)
            % Computes a smooth spline of a single graph edge.
            % Based on original function by Kin Quan 2018
            
            % The input is the list of ordered index
            % The output is the smooth spline as a matlab sturct
            
            % get linear indexed points of current and previous branch,
            % combine.
            previous_awy = find([obj.Glink(:).n2] == obj.Glink(link_index).n1);
            [x_p1, y_p1, z_p1] = ind2sub(size(obj.CT), obj.Glink(previous_awy).point);
            [x_p2, y_p2, z_p2] = ind2sub(size(obj.CT), obj.Glink(link_index).point);
            x_point = [x_p1, x_p2];
            y_point = [y_p1, y_p2];
            z_point = [z_p1, z_p2];

            %Smooth all points using moving average
            voxel_sz = obj.CTinfo.PixelDimensions;
            smooth_x = smooth(x_point*voxel_sz(1),11, 'moving');
            smooth_y = smooth(y_point*voxel_sz(2),11, 'moving');
            smooth_z = smooth(z_point*voxel_sz(3),11, 'moving');
            
            % extract just current airway smoothed points
            csmooth_x = smooth_x(length(x_p1)+1:end);
            csmooth_y = smooth_y(length(x_p1)+1:end);
            csmooth_z = smooth_z(length(x_p1)+1:end);
            
            %Complete smooth data
            smooth_data_points = [csmooth_x csmooth_y csmooth_z]';
            
            %Generating the spline
            obj.Splines{link_index, 1} = cscvn(smooth_data_points);
                end
        
        function t_points = ComputeSplinePoints(obj, link_index)
            % * Compute Spline if necessary
            if isempty(obj.Splines{link_index, 1})
                ComputeSpline(obj, link_index)
            end
            spline = obj.Splines{link_index, 1};
            
            [t_points, arc_length] = spline_points(spline, obj.spline_sample_sz);
            
            % save parametrised sample points in second column.
            obj.Splines{link_index, 2} = t_points;
            
            % save arc_length measurement of each sample point.
            obj.arclength{link_index, 1} = arc_length;
            
    end
        
        %% TRAVERSING AIRWAYS METHODS %%%
        % Technical: These methods traverse the airway centreline spline 
        % interpolating the CT image at each spline sample point.
        function obj = CreateAirwayImage(obj, link_index)
            % Constructs perpendicular images as if travelling along an
            % airway segment in CT image and Segmentation.
            spline_t_points = ComputeSplinePoints(obj, link_index);
            % set up slice store
            TransAirwayImage = cell(length(spline_t_points),1);
            TransSegImage = cell(length(spline_t_points),1);
            for i = 1:length(spline_t_points)
                spline = obj.Splines{link_index, 1};
                % * Compute Normal Vector per spline point
                [normal, CT_point] = AirQuant.ComputeNormal(spline, ...
                    spline_t_points(i));
                % get approx airway size from distance map
                approx_diameter = ComputeDmapD(obj, CT_point);
                % compute intepolated slice size
                plane_sz = ceil(approx_diameter*obj.plane_scaling_sz);
                % use max plane size if current plane size exceeds it
                if plane_sz > obj.max_plane_sz
                    plane_sz = obj.max_plane_sz;
                end
                % * Interpolate Perpendicular Slice per spline point
                [InterpAirwayImage, InterpSegImage] = ...
                    InterpolateCT(obj, normal, CT_point,  ...
                    plane_sz, obj.plane_sample_sz);
                % * Replace NaN entries in images with zero.
                InterpAirwayImage(isnan(InterpAirwayImage)) = 0;
                InterpSegImage(isnan(InterpSegImage)) = 0;
                TransAirwayImage{i,1} = InterpAirwayImage;
                TransSegImage{i,1} = InterpSegImage;
            end
            % * Save traversed image and arclength for each image
            obj.TraversedImage{link_index, 1} = TransAirwayImage;
            obj.TraversedSeg{link_index, 1} = TransSegImage;
            % save obj to disk after every branch
            save(obj)
        end
        
        
        function [CT_plane, seg_plane] = InterpolateCT(obj, normal, ...
                CT_point,  plane_sz, plane_sample_sz)
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
                basis_vecs(:,2), CT_point, plane_sz, plane_sample_sz);
            
            % * Execute cubic inperpolation on CT
            plane_intensities = interp3(x_domain,y_domain,z_domain,...
                obj.CT,plane_grid.y(:),plane_grid.x(:),...
                plane_grid.z(:),'cubic');
            
            % * Execute cubic inperpolation on seg
            seg_intensities = interp3(x_domain,y_domain,z_domain,...
                double(obj.seg),plane_grid.y(:),plane_grid.x(:),...
                plane_grid.z(:),'cubic');
            
            % Reshape
            plane_length = sqrt(length(plane_grid.y(:)));
            CT_plane = reshape(plane_intensities,...
                [plane_length plane_length]);
            seg_plane = reshape(seg_intensities,...
                [plane_length plane_length]);
        end
        
        function approx_diameter = ComputeDmapD(obj, CT_point)
            % Convert CT_point mm back to voxel ind
            vox_point = CT_point'./obj.CTinfo.PixelDimensions;
            % find nearest skeleton point to voxpoint
            if isempty(obj.skel_points)
                ComputeSkelPoints(obj)
            end
            P = obj.skel_points;
            k = dsearchn(P,vox_point);
            % Get radius and convert to diameter
            approx_diameter = obj.Dmap(P(k,1),P(k,2),P(k,3))*2;
            % incase of edge case, unit radius
            if approx_diameter == 0
                approx_diameter = 2;
            end
            % TODO: diameter conversion to mm
        end
        
        %% MEASUREMENT METHODS
        % Methods that measure the airway size on interpolated CT images.
        function obj = FindAirwayBoundariesFWHM(obj, link_index)
            %Based on function by Kin Quan 2018 that is based on Kiraly06
            
            slices_sz = size(obj.TraversedImage{link_index, 1}, 1);
            
            % Prepping the outputs
            raycast_FWHMl = cell(slices_sz, 1);
            raycast_FWHMp = cell(slices_sz, 1);
            raycast_FWHMr = cell(slices_sz, 1);
            
            % For every traversed slice
            for k = 1:slices_sz
                try
                    % * Compute airway centre
                    % Check that airway centre is slice centre
                    center = ...
                        Return_centre_pt_image(...
                        obj.TraversedSeg{link_index, 1}{k,1});
                    
                    % Recompute new centre if necessary
                    [centre_ind , new_centre] =  ...
                        Check_centre_with_segmentation(...
                        obj.TraversedImage{link_index, 1}{k,1}, ...
                        obj.TraversedSeg{link_index, 1}{k,1});
                    if ~centre_ind
                        center = fliplr(new_centre);
                    end
                    
                    % * Raycast
                    [CT_rays, seg_rays, coords]= Raycast(obj, ...
                        obj.TraversedImage{link_index, 1}{k,1}, ...
                        obj.TraversedSeg{link_index, 1}{k,1}, center);
                    
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
            image_sz = size(interpslice);
            limit_row = abs(image_sz(1) - center(1));
            limit_col = abs(image_sz(2) - center(2));
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
                pi*obj.plane_sample_sz.^2;
        end
        
        %% EXTENT Methods
        % Analysis: Methods looking at the extent of airway existance that
        % the segmentation covers.
        function counts = searchgenperlobe(obj, gen, lobe)
            Nawy = length(obj.Glink);
            counts = 0;
            for i = 1:Nawy
                if obj.Glink(i).generation == gen && strcmp(obj.Glink(i).lobe, lobe)
                counts = counts + 1;
                end
            end
        end
        
        function genperlobe = extentstats(obj)
            % make table of lobes
            maxgen = max([obj.Glink(:).generation]);
            sz = [maxgen, 6];
            varTypes = cell(1, 6);
            varTypes(:) = {'double'};
            varNames = AirQuant.LobeLabels();
            genperlobe = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
            % get number of airways of that generation in that lobe
            for lobe = 1:6
                lobeid = varNames{1,lobe};
                for gen = 1:maxgen
                    genperlobe.(lobeid)(gen) = searchgenperlobe(obj, gen, lobeid);
                end
            end
        end
        
        function [overextentidx, compare] = compareextent(obj, populationstats, threshold)
            % threshold = minimum gen expected
            % populationstats = table generated by summarising results from
            % genperlobe over several images.
            subject_genperlobe = table2array(extentstats(obj));
            if istable(populationstats)
                populationstats = table2array(populationstats);
            end
            gendiff = size(populationstats, 1) - size(subject_genperlobe, 1);
            if gendiff < 0
                populationstats = [populationstats; zeros(abs(gendiff),6)];
            elseif gendiff > 0
                subject_genperlobe = [subject_genperlobe; zeros(abs(gendiff),6)];
            end
                
            % compare this subject to population
            compare = floor(subject_genperlobe - populationstats);
            
            overextentidx = false(length(obj.Glink),1);
%             safenodes = [];
            varNames = AirQuant.lobelabels();
            % BF search
            [~,bfs_Gind] = bfsearch(obj.Gdigraph, 1,'edgetonew');
            % convert graph idx into link index
            bfs_Lind = zeros(size(bfs_Gind));
            for ind = 1:length(bfs_Gind)
                bfs_Lind(ind) = obj.Gdigraph.Edges.Label(bfs_Gind(ind));
            end
            % reverse order
            bfs_Lind = flip(bfs_Lind);
            for lobe = 1:6
                lobeid = varNames{1,lobe};
                for gen = size(subject_genperlobe, 1):-1:threshold+1
                    n_overextent = compare(gen, lobe);
                    if n_overextent <= 0
                        continue
                    end
                    overextentawy(obj, gen, lobeid, n_overextent);
                end
            end
%             for awy = 1:length(obj.Glink)
%                 if ~isempty(ismember(obj.Glink(awy).generation,[1:threshold]))
%                     overextentidx(awy) = false;
%                 end
%             end
%             function overextentawy(obj, gen, lobe, n_overextent)
%             Nawy = length(obj.Glink);
%             counts = 1;
%             mode = 0;
%             for i = 1:Nawy
%                 if counts > n_overextent
%                     mode = 1; % safe mode
%                 end
%                 checksafe = isempty(ismember([obj.Glink(i).n1, obj.Glink(i).n2], safenodes));
%                 if obj.Glink(i).generation == gen && strcmp(obj.Glink(i).lobe, lobe) && checksafe == 0
%                     if mode == 0
%                         overextentidx(i) = true; % highlight as excess
%                         counts = counts + 1;
%                     elseif mode == 1
%                         safenodes = [safenodes; obj.Glink(i).n1; obj.Glink(i).n2];
%                     end
%                 end
%             end
%             end    
%         end

            function overextentawy(obj, gen, lobe, n_overextent)
            counts = 1;
            for i = 1:length(bfs_Lind)
                if counts > n_overextent
                    break
                end
                awyind = bfs_Lind(i);
                if obj.Glink(awyind).generation == gen && strcmp(obj.Glink(awyind).lobe, lobe)
                   overextentidx(awyind) = true; % highlight as excess
                   counts = counts + 1;
                end
            end
            end    
        end
         
        function h = plotoverextent(obj, overextentidx, type)
            overextentedge = overextentidx(obj.Gdigraph.Edges.Label);
            h = plot(obj, type);
            highlight(h,'Edges',overextentedge, 'LineStyle', ':', 'LineWidth',1)
        end
     
        %% TAPERING ANALYSIS METHODS
        % Analysis: Airway tapering metrics.
        
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
            path = shortestpath(G, obj.carina_node, terminal_node_idx);
            
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
                cum_arclength = [cum_arclength; current_arclength_array'];
                
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
            terminalnodelist(terminalnodelist ==  1) = [];
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
%             if nargout > 1
            obj.Specs.AllTaperResults = AllTaperResults;
%             end
        end
        
        function [intrataper, averagediameter] = ComputeIntraTaperAll(obj, prunelength)
            if nargin < 2
                prunelength = [0 0];
            end
            % loop through branches
            intrataper = NaN(length(obj.arclength), 3);
            averagediameter = NaN(length(obj.arclength), 3);
            for ii = 1:length(obj.arclength)
                if ii == obj.trachea_path
                    continue
                end
                
                [intrataper(ii,:), averagediameter(ii,:)] = ComputeIntraTaper(obj, prunelength, ii);
            end
        end
        
        function [intrataper, averagediameter] = ComputeIntraTaper(obj, prunelength, idx, plotflag)
            % prunelength given as a 2 element vector, the length in mm to
            % ignore at the begining and end of the branch.
            
            if nargin < 4
                plotflag = 0;
            end
            
            intrataper = NaN(1, 3);
            averagediameter = NaN(1, 3);
            
            % get arclength
            al = obj.arclength{idx, 1};
            
            % get branch radii
            areas = zeros(length(obj.FWHMesl{idx,1}), 3);
            for jj = 1:length(areas)
                try % incase area is NaN
                    areas(jj,1) = obj.FWHMesl{idx, 1}{jj, 1}.area;
                    areas(jj,2) = obj.FWHMesl{idx, 2}{jj, 1}.area;
                    areas(jj,3) = obj.FWHMesl{idx, 3}{jj, 1}.area;
                catch
                    areas(jj,1) = NaN;
                    areas(jj,2) = NaN;
                    areas(jj,3) = NaN;
                end
            end
            
            % prune ends
            disp(idx)
            prune = (al >= prunelength(1) & al <= al(end) - prunelength(2));
            al = al(prune);
%             areas = areas(repmat(prune',1,3));
            coeff = NaN(2,3);
            % convert area to diameters
            diameters = sqrt(areas/pi)*2;
            for jj = 1:3
                try % incase no branch left after pruning/too few points
                    Dvec = diameters(prune,jj);
                    % fit bisquare method
                    coeff(:,jj) = robustfit(al, Dvec,'bisquare');
                    % compute intra-branch tapering as percentage
                    intrataper(jj) = -coeff(2,jj)/coeff(1,jj) * 100;
                    % compute average area
                    averagediameter(jj) = mean(Dvec, 'omitnan');
                catch
                    % leave as NaN
                end
            end
            
            if plotflag == 1 && ~any(isnan(intrataper))
                titlevec = ["inner"; "peak"; "outer"];
                for jj = 1:3
                subplot(3,1,jj)
                plot(al, areas(prune,jj), 'k.')
                hold on
                plot(al, coeff(2,jj)*al + coeff(1,jj),'r')
                legend('data', 'bisquare fit', 'Location', 'best')
                xlabel('arclength (mm)')
                ylabel('area mm^2')
                title(sprintf('Branch idx: %i %s intrataper value: %.2f%% average: %.2f mm.',...
                idx, titlevec(jj), intrataper(jj), averagediameter(jj)))
                hold off
                
                end
            end
            
            end
        
        function intertaper = ComputeInterTaper(obj, prunelength)
            if nargin < 2
                prunelength = [0 0];
            end
            % use output from intrataperall
            [~, averagediameter] = ComputeIntraTaperAll(obj, prunelength);
            % loop through branches
            intertaper = NaN(length(obj.arclength), 3);
            for ii = 1:length(averagediameter)
                if ii == obj.trachea_path
                    continue
                end
                for jj = 1:3
                    % identify parent by predecessor node
                    parent = find([obj.Glink.n2] == obj.Glink(ii).n1);
                    intertaper(ii,jj) = (averagediameter(parent, jj) - averagediameter(ii,jj))...
                        /(averagediameter(parent, jj)) * 100;
                end
            end

        end
        
        function SegmentTaperResults = SegmentTaperAll(obj, prunelength)
            
            % compute taper results by segment
            [intrataper, avg] = ComputeIntraTaperAll(obj, prunelength);
            intertaper = ComputeInterTaper(obj, avg);
            
            % organise into column headings
            branch = 1:length(obj.Glink);
            
            inner_intra = intrataper(:, 1);
            peak_intra = intrataper(:, 2);
            outer_intra = intrataper(:, 3);
            
            inner_avg = avg(:, 1);
            peak_avg = avg(:, 2);
            outer_avg = avg(:, 3);
            
            inner_inter = intertaper(:, 1);
            peak_inter = intertaper(:, 2);
            outer_inter = intertaper(:, 3);
            
            % convert to table
            SegmentTaperResults = table(branch', inner_intra, peak_intra, ...
                outer_intra, inner_avg, peak_avg, outer_avg, ...
                inner_inter, peak_inter, outer_inter);
            
            % add gen info
            SegmentTaperResults.generation = [obj.Glink.generation]';
            
            % add lobe info if available
            try % only add lobe information if it exists
                SegmentTaperResults.lobe = [obj.Glink.lobe]';
                % sort by lobe
                SegmentTaperResults = sortrows(SegmentTaperResults, 'lobe');
            catch
            end
            
            % Save to AQ object
            obj.Specs.SegmentTaperResults = SegmentTaperResults;
        end
        
        %% TAPERING VISUALISATION METHODS
        % Visualisation: viewing results of taper metrics.
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
        
        function TaperBoxPlot(obj, type)
            % requires lobe classification
            if ~isfield(obj.Glink, 'lobe')
                error('Lobe classification required. Please run ComputeAirwayLobes() first.')
            end
            % run ComputeTaperAll if analysis not saved.
            if  ~isfield(obj.Specs, 'AllTaperResults')
                AllTaperResults = ComputeTaperAll(obj);
            else
                AllTaperResults = obj.Specs.AllTaperResults;
            end   
            if nargin < 2
                type = 'other';
            end
                        
            switch type
                case 'inner'
                    logtaperdata = [AllTaperResults.logtapergrad_inner];
                    typetxt = 'Inner lumen';
                case 'peak'
                    logtaperdata = [AllTaperResults.logtapergrad_peak];
                    typetxt = 'Peak wall';
                case 'outer'
                    logtaperdata = [AllTaperResults.logtapergrad_outer];
                    typetxt = 'Outer wall';
                otherwise
                    subplot(3,1,1)
                    TaperBoxPlot(obj, 'inner')
                    subplot(3,1,2)
                    TaperBoxPlot(obj, 'peak')
                    subplot(3,1,3)
                    TaperBoxPlot(obj, 'outer')
                    return
            end
            datalabels = [AllTaperResults.lobe];
            
            boxplot(logtaperdata, datalabels)
            xlabel('Lobe')
            ylabel('Log Taper Gradient')
            title(typetxt)
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
        
        %% VISUALISATION
        %%% Airway Strucutral Tree
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
                case {'lobe','lobes'}
                    try
                        lobes = [obj.Glink(:).lobe];
                        edgelabels = lobes(G.Edges.Label);
                    catch
                        error('Need to run ComputeAirwayLobes first')
                    end
                case {'generation','gen'}
                    gens = [obj.Glink(:).generation];
                    edgelabels = gens(G.Edges.Label);
                case 'none'
                    edgelabels = '';
                otherwise
                    warning('Unexpected plot type. Resorting to default type.')
                    edgelabels = G.Edges.Label;
            end
            
            h = plot(G,'EdgeLabel',edgelabels, 'Layout', 'layered');
%             h = plot(G, 'Layout', 'layered');
            h.NodeColor = 'r';
            h.EdgeColor = 'k';
            % check for error edges
            [~,~,erroredge] = DebugGraph(obj);
            if ~isempty(erroredge)
                highlight(h,'Edges',erroredge,'EdgeColor','r');
            end
        end

        function PlotTree(obj, gen)
            % Plot the airway tree with nodes and links in image space. Set
            % gen to the maximum number of generations to show.
            % Original Function by Ashkan Pakzad on 27th July 2019.
                            
            if nargin == 1
                gen = max([obj.Glink(:).generation]);
            end
            
%             isosurface(obj.skel);
%             alpha(0.7)
            % set up reduced link graph and skel 
            vis_Glink_logical = [obj.Glink(:).generation] <= gen;
            vis_Glink_ind = find(vis_Glink_logical == 1);
            vis_Glink = obj.Glink(vis_Glink_logical);
            % set up reduced node graph
            vis_Gnode_logical = false(length(obj.Gnode),1);
            for i = 1:length(obj.Gnode)
                if ~isempty(intersect(obj.Gnode(i).links, vis_Glink_ind))
                    vis_Gnode_logical(i) = 1;
                end
            end
            vis_Gnode_ind = find(vis_Gnode_logical == 1);
            vis_Gnode = obj.Gnode(vis_Gnode_logical);
            % set up reduced skel
            vis_skel = zeros(size(obj.skel));
            vis_skel([vis_Glink.point]) = 1;
            %vis_skel([vis_Gnode.idx]) = 1;
            
            isosurface(vis_skel)
            alpha(0.7)
            
            hold on
            
            % edges
            ind = zeros(length(vis_Glink), 1);
            for i = 1:length(vis_Glink)
                ind(i) = vis_Glink(i).point(ceil(end/2));
            end
            [Y, X, Z] = ind2sub(size(obj.skel),ind);
            nums_link = string(vis_Glink_ind);
            %plot3(X,Y,Z, 'b.', 'MarkerFaceColor', 'none');
            text(X+1,Y+1,Z+1, nums_link, 'Color', [0, 0.3, 0])
            
            % nodes
            X_node = [vis_Gnode.comy];
            Y_node = [vis_Gnode.comx]; 
            Z_node = [vis_Gnode.comz];
            nums_node = string(vis_Gnode_ind);
            plot3(X_node,Y_node,Z_node, 'r.', 'MarkerSize', 18, 'Color', 'r');
            text(X_node+1,Y_node+1,Z_node+1, nums_node, 'Color', [0.8, 0, 0])
            
            %axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            view(80,0)
            axis vis3d
            % undo matlab display flip
            ax = gca;
            ax.XDir = 'reverse';
            
        end
        
        %%% Splines
        function PlotSplineTree(obj)
            % loop through every branch, check spline has already been
            % computed, compute if necessary. Skip trachea. Plot spline.
            % TODO: This may not work as well when VoxelSize ~= [1,1,1]
            for i = 1:length(obj.Glink)
                if isempty(obj.Splines{i, 1})
                    if ismember(i, obj.trachea_path)
                        continue
                    else
                        ComputeSpline(obj, i)
                    end
                end
                fnplt(obj.Splines{i, 1})
                hold on
            end
        end
        
        function h = PlotSplineVecs(obj, subsamp, link_index)
            % Plot and Visualise tangental vectors off spline sample
            % points. 1/subsamp = proportion of spline points to sample,
            % default subsamp = 2. link_index = airway indices to plot,
            % default link_index = all.
            
            % plot all branches if specific airway not provided
            if nargin < 2
                subsamp = 2;
            end
            
            if nargin < 3
                link_index = 1:length(obj.Glink);
            end
            
            for iidx = link_index
                % generate spline and points if it doesnt exist
                if isempty(obj.Splines{iidx, 1})
                    ComputeSplinePoints(obj, iidx);
                end
                
                % get vecs and origin for spline
                spline = obj.Splines{iidx, 1};
                samplepnts = obj.Splines{iidx, 2};
                vecs = zeros(3, length(samplepnts));
                origins = zeros(3, length(samplepnts));
                for jj = 1:length(samplepnts)
                    point = samplepnts(jj);
                    [vecs(:,jj), origins(:,jj)]= AirQuant.ComputeNormal(spline, point);
                end
                
                % subsample data, i.e. delete a portion
                vecs = vecs(:,1:subsamp:end); origins = origins(:,1:subsamp:end);
                
                % rescale origins from mm to vox and swap x and y in plot.
                origins = origins./obj.CTinfo.PixelDimensions';
                vecs = vecs./obj.CTinfo.PixelDimensions';

                
                % plot vectors with translucent airway volume
                h = quiver3(origins(2,:),origins(1,:),origins(3,:),...
                    vecs(2,:),vecs(1,:),vecs(3,:));
                branch_seg = obj.ClassifySegmentation(); % get labelled segmentation
                patch(isosurface(branch_seg == iidx), ...
                    'FaceAlpha', 0.3, 'FaceColor', [0 .55 .55], 'EdgeAlpha', 0);
                hold on
            end
            
            % visualise the connecting airway segs if only 1 airway
            % requested
            if length(link_index) < 2
                % find connecting airways
                nodes = [obj.Glink(link_index).n1, obj.Glink(link_index).n2];
                conn_awy = [];
                for connii = 1:length(obj.Glink)
                    for nodeii = nodes
                        conn_awy = [conn_awy find([obj.Glink.n1] == nodeii | ...
                            [obj.Glink.n2] == nodeii)];
                    end
                end
                conn_awy = unique(conn_awy);
                conn_awy(conn_awy == link_index) = [];
                % visualise them
                for awyii = conn_awy
                    patch(isosurface(branch_seg == awyii), ...
                        'FaceAlpha', 0.3, 'FaceColor', [.93 .79 0], 'EdgeAlpha', 0);
                end
            end
            vol3daxes(obj)
            hold off
        end
        
        %%% Volumetric
        function PlotMap3D(obj, mode)
            % Recommend to use View3D if colour labels appear buggy.
            
            % mode = 'TaperGradient', 'generation', 'lobes'
            axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            
            % generating the color data
            cdata = zeros(size(obj.seg));
            branch_seg = ClassifySegmentation(obj);
            switch mode 
                case 'tapergradient'
                    % TODO: rewrite this bit....
                    for i = 1:length(obj.Specs)
                        cdata(branch_seg == i) = obj.Specs(i).FWHMl_logtaper*-1;
                    end
                case 'generation'
                    for i = 1:length(obj.Glink)
                        % add 1 to gen index to differentiate from bg 0.
                        cdata(branch_seg == i) = obj.Glink(i).generation+1;
                    end
                    clims = [1 max(cdata(:))];
                    colourshow = clims(1):clims(2);
                    colorbarstring = 'Generation Number';
                    % reduce colourlabels by 1 from cdata to reflect true
                    % gen.
                    colourlabels = 0:max(cdata(:))-1;
                    maptype = 'sequential';

                case 'lobe'
                    % convert lobe id to number
                    lobeid = {'B','RU','RM','RL','LU','LUlin','LL'};
                    for i = 1:length(obj.Glink)
                        cdata(branch_seg == i) = find(strcmp(lobeid, obj.Glink(i).lobe));
                    end
                    clims = [1 max(cdata(:))];
                    colorbarstring = 'Lobe';
                    colourshow = clims(1):clims(2);
                    colourlabels = lobeid;
                    maptype = 'qualitative';
                    
                otherwise
                    error('Choose appropiate mode.')
            end
            % produce segmentation 3d object
            p = patch(isosurface(obj.seg));
            
            %%% assign vertex face colour index by nearest point on volume.
            % get volume points
            cdata_pnt = find(cdata > 0);
            % convert to list of subindices
            [y,x,z] = ind2sub(size(cdata), cdata_pnt);
            % search for nearest point of each vertex origin
            near_i = dsearchn([x,y,z], p.Vertices);
            % assign colour index to that vertex
            p.FaceVertexCData = cdata(cdata_pnt(near_i));
            p.FaceColor = 'flat';
            p.EdgeColor = 'none';
            
            % set up colourmap
            map = linspecer(max(cdata(:)), maptype);
            colormap(map)
            c = colorbar('Ticks', colourshow, 'TickLabels', colourlabels);
            c.Label.String = colorbarstring;
            caxis(clims)
            vol3daxes(obj)
        end
        
        function View3D(obj, mode)
            % View segmentation volume with different labels. In MATLAB's
            % volviewer.
            % mode = 'TaperGradient', 'generation', 'lobes'
            
            % generating the color data
            labelvol = zeros(size(obj.seg));
            branch_seg = ClassifySegmentation(obj);
            switch mode 
                case 'generation'
                    for i = 1:length(obj.Glink)
                        labelvol(branch_seg == i) = obj.Glink(i).generation;
                    end
                case 'lobe'
                    % convert lobe id to number
                    lobeid = {'B','RU','RM','RL','LU','LUlin','LL'};
                    for i = 1:length(obj.Glink)
                        labelvol(branch_seg == i) = find(strcmp(lobeid, obj.Glink(i).lobe))-1;
                    end
                otherwise
                    error('Choose appropiate mode. "Generation" or "Lobe".')
            end
            
            % undo matlabs X-axis flip for viewing.
            labelvol = flip(labelvol,1);
            seg_view = flip(obj.seg,1);
            
            % Generate suitable label colours
            map = linspecer(max(labelvol(:))+1);

            % vol viewer with labels to display.
            figure;
            labelvolshow(labelvol, seg_view, ...
                'ScaleFactors', obj.CTinfo.PixelDimensions, ...
                'LabelColor', map, 'BackgroundColor', [1,1,1], ...
                'CameraPosition', [-4.2 0.8  2], 'CameraViewAngle', 10, ...
                'CameraTarget', [0, 0, 0.1]);
        end
        
        function PlotSegSkel(obj)
        % plot segmentation and skeleton within each other.
        patch(isosurface(obj.seg),'EdgeColor', 'none','FaceAlpha',0.3);
        hold on
        isosurface(obj.skel)
        vol3daxes(obj)
        
        end
        
        function vol3daxes(obj, ax)
            % utility function for 3D volumetric plotting. Sets the aspect
            % ratio according to voxel size and reverses the x axes for LPS
            % viewing.
            
            if nargin < 2 % current axes if not specified
                ax = gca;
            end
            axis vis3d
            view(80,0)
            % aspect ratio
            ax.DataAspectRatio = 1./obj.CTinfo.PixelDimensions;
            % undo matlab display flip
            ax.XDir = 'reverse';
        end
        
        %%% CT Airway slices
        function PlotAirway3(obj, link_index)
            % Plot resampled airway slices overlayed with FWHMesl ray cast
            % points and fitted ellipse
            f = figure('Position',  [100, 100, 850, 600]);
            slide = 1;
            PlotAirway(obj, link_index, slide)
            numSteps = size(obj.TraversedImage{link_index,1}, 1);
            
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
            canvas_sz = floor(obj.max_plane_sz/obj.plane_sample_sz);
            canvas = nan(canvas_sz);
            image = obj.TraversedImage{link_index, 1}{slide,1};
            image_sz = size(image,1);
            min_centre = canvas_sz/2 - image_sz/2;
            max_centre = canvas_sz/2 + image_sz/2;
            canvas(min_centre+1:max_centre, min_centre+1:max_centre) = image;
            imagesc(canvas)
            colormap gray
            
            try % try block incase FWHMesl has not been executed.
                % plot ray cast results
                
%                 plot(obj.FWHMesl{link_index, 1}{slide, 1}.x_points, obj.FWHMesl{link_index, 1}{slide, 1}.y_points,'r.')
%                 plot(obj.FWHMesl{link_index, 2}{slide, 1}.x_points, obj.FWHMesl{link_index, 2}{slide, 1}.y_points,'c.')
%                 plot(obj.FWHMesl{link_index, 3}{slide, 1}.x_points, obj.FWHMesl{link_index, 3}{slide, 1}.y_points,'y.')
                
                % plot ellipse fitting
                ellipse(obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(4),...
                    obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(1)+min_centre,...
                    obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(2)+min_centre,'m');
                
%                 ellipse(obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(4),...
%                     obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(1),...
%                     obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(2),'b');
                
                ellipse(obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(4),...
                    obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(1)+min_centre,...
                    obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(2)+min_centre,'y');
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
                text(ax, 1,5,sprintf('Arc Length = %4.2f mm; Inner area = %4.2f mm^2; Peak area = %4.2f mm^2; Outer area = %4.2f mm^2; %3.0i of %3.0i', ...
                    obj.arclength{link_index, 1}(slide), obj.FWHMesl{link_index, 1}{slide, 1}.area, obj.FWHMesl{link_index, 2}{slide, 1}.area ,...
                    obj.FWHMesl{link_index, 3}{slide, 1}.area, slide, size(obj.TraversedImage{link_index, 1},3)));
            catch
                text(ax, 1,5,sprintf('Arc Length = %4.1f mm; %3.0i of %3.0i', ...
                    obj.arclength{link_index, 1}(slide), slide, size(obj.TraversedImage{link_index, 1},3)));
            end
        end
        
        function s = OrthoViewAirway(obj, link_index)
        % View a series of an airway segment's slices as a volume image 
        % stack using MATLAB's inbuilt othogonal 3d viewer.
        
        % convert from cell stack to 3D array.
        awycell =  obj.TraversedImage{link_index,1};
        canvas_sz = floor(obj.max_plane_sz/obj.plane_sample_sz);
        awyarray = zeros([canvas_sz, canvas_sz, length(awycell)]);
        for slice = 1:length(awycell)
            image = obj.TraversedImage{link_index, 1}{slice,1};
            image_sz = size(image,1);
            min_centre = canvas_sz/2 - image_sz/2;
            max_centre = canvas_sz/2 + image_sz/2;
            awyarray(min_centre+1:max_centre, min_centre+1:max_centre, slice) = image;
        end
        % display with orthoview
        fig = figure;
        s = orthosliceViewer(awyarray, 'DisplayRangeInteraction','off', ...
            'ScaleFactors',[obj.plane_sample_sz, obj.plane_sample_sz, obj.spline_sample_sz],...
            'CrosshairLineWidth', 0.3);
        % Can only alter size of figure window after orthosliceviewer.
        fig.Name = ['AirQuant: Airway Ortho View. Idx ', mat2str(link_index), '.'];
        fig.Units = 'normalized';
        fig.Position = [0.1,0.01,0.6,0.9];
        end
        
        %% EXPORT METHODS
        % methods for exporting processed data.
        function exportlobes(obj, savename)
            % export airway segmentation labelled by lobes to nii.gz
            
            % classify by branch first
            branch_seg = ClassifySegmentation(obj);
            lobeid = {'B','RU','RM','RL','LU','LUlin','LL'};
            lobe_airway_seg = zeros(size(obj.seg));
            
            % convert branch classification to lobe classification
            try
                for i = 1:length(obj.Glink)
                    lobe_airway_seg(branch_seg == i) = find(strcmp(lobeid, obj.Glink(i).lobe));
                end
            catch
                error('Need to run ComputeAirwayLobes first')
            end
            
            % reduce datatype and change header info
            lobe_airway_seg = uint8(lobe_airway_seg);
            header = obj.CTinfo;
            header.Datatype = 'uint8';
            header.BitsPerPixel = '8';
            header.Description = 'Airway Lobe Segmentation using AirQuant, layers: B,RU,RM,RL,LU,LUlin,LL';
            
            niftiwrite(lobe_airway_seg, savename, header, 'Compressed', true);
        end
        
        function innerareas = GetAreas(obj)
            % Extract the inner area of all branches from FWHM results and output
            
            if  ~isfield(obj.Specs, 'innerareas')
                innerareas = cell(length(obj.Glink), 1);
                for ii = 1:length(obj.Glink)
                    for jj = 1:length(obj.FWHMesl{ii,1})
                        try
                            innerareas{ii, 1} = [innerareas{ii, 1}, obj.FWHMesl{ii,1}{jj, 1}.area];
                        catch
                            innerareas{ii, 1} = [innerareas{ii, 1}, NaN];
                        end
                    end
                end
                maxlen = max(cellfun(@length, innerareas));
                innerareas = cellfun(@(x)([x nan(1, maxlen - length(x))]), innerareas, 'UniformOutput', false);
                innerareas = cell2mat(innerareas);
                obj.Specs.innerareas = innerareas;
            else
                innerareas = obj.Specs.innerareas;
            end
            
        end
        
        function innerdiameters = GetDiameters(obj)
            % Extract the inner area of all branches from FWHM results,
            % convert to diameter and output.
            if  ~isfield(obj.Specs, 'innerdiameters')
                innerareas = GetAreas(obj);
                innerdiameters = sqrt(innerareas/pi)*2;
                obj.Specs.innerdiameters = innerdiameters;
            else
                innerdiameters = obj.Specs.innerdiameters;
            end
        end
        
        function exportraw(obj, savename)
            % export Glink table and arclengths and inner diameter
            % measurements
            
            % Glink to table
            BranchInfo = struct2table(obj.Glink);
            
            % arclengths to table
            arclengths = obj.arclength;
            maxlen = max(cellfun(@length, arclengths));
            arclengths = cellfun(@(x)([x nan(1, maxlen - length(x))]), arclengths, 'UniformOutput', false);
            arclengths = cell2mat(arclengths);
            
            % diameters
            innerdiameters = GetDiameters(obj);
            
            % export csvs
            writetable(BranchInfo, ['BranchInfo_',savename, '.csv']);
            writematrix(arclengths, ['ArcLengths_',savename, '.csv']);
            writematrix(innerdiameters, ['InnerDiameters_',savename, '.csv']);

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
        
        function labels = LobeLabels()
            % returns cell of lobe lables used in AirQuant
            labels = {'RU','RM','RL','LU','LUlin','LL'};
        end
    end
end