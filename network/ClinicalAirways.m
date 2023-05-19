classdef ClinicalAirways < TubeNetwork
    % short desc
    %
    % long desc
    %
    % .. todo: add documentation to this function
    %
    % Args:
    %   x(type):
    %
    % Return:
    %   y(type):
    %

    methods
        function obj = ClinicalAirways(varargin)
            % short desc
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            obj@TubeNetwork(varargin{:});
            
            obj.regioncategories.lobe = {'B','RUL','RML','RLL','LUL','LML','LLL','T'};
            obj.IdentifyCarinaAndTrachea();
            
            % attempt to classify into lunglobes
            try
                obj.ClassifyLungLobes()
            catch
            end

        end

        function obj = MakeTubes(obj, glink)
            for ii = 1:length(glink)
                obj.tubes = [obj.tubes, Airway(obj, glink(ii).point, ii)];
            end
        end

        function obj = IdentifyCarinaAndTrachea(obj)
            % short desc
            %
            % long desc
            %
            % .. todo:
            %   * May need to add further methods to reidentify generations
            %   for all airways.
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %

            % identfy carina
            g = TubesAsEdges(obj);
            [~, carina_node] = max(centrality(g,'outcloseness'));
            carinaend_gidx = inedges(g,carina_node);
            obj.tubes(g.Edges.ID(carinaend_gidx)).SetCarinaEnd();

            % all tubes predecessing the carina is considered part of the
            % trachea and are therefore generation 0
            tracheanodes = predecessors(g, carina_node);
            for nid = tracheanodes
                eid = outedges(g, nid);
                tubeid = g.Edges.ID(eid);
                for tid = tubeid'
                    obj.tubes(tid).SetTrachea();
                end
            end

            % reclass all tube generations by n decendants from 0 gen.
            obj.RunAllTubes('SetGeneration');
        end

        function obj = ClassifyLungLobes(obj)
            % Attempt to classify airways into their lobes.
            %
            % Adapted from original algorithm by Gu et al. 2012.
            %
            %

            % find trachea tubes and assign
            tubegens = [obj.tubes(:).generation];
            tracheaID = find(tubegens == 0);
            for ii = tracheaID
                obj.tubes(ii).SetRegion('lobe','T');
                obj.tubes(ii).SetRegion('name','Trachea');
            end

            % find carina-end tube and init algorithm
            carinaendID = find([obj.tubes(:).carinaend]);
            assert(length(carinaendID) == 1, ['There should only be one' ...
                ' carina end. This indicates a major failure in the lobe ' ...
                'classification algorithm.'])

            % Identify left and right major bronchi by comparing ends of
            % child of carina-end trachea
            MB = obj.tubes(carinaendID).children;
            [mX, ~, ~] = obj.I2S(ClinicalAirways.SkelEnds(MB));
            [~,leftMBI] = max(mX);
            [~,rightMBI] = min(mX);
            MB(leftMBI).SetRegion('lobe','B');
            MB(leftMBI).SetRegion('name','LeftMajor');
            MB(rightMBI).SetRegion('lobe','B');
            MB(rightMBI).SetRegion('name','RightMajor');

            % Identify upper/lingular and lower left lobe 'LLL'
            MLlung = MB(leftMBI).children;
            [~, ~, MLlz] = obj.I2S(ClinicalAirways.SkelEnds(MLlung));
            [~,MLULLML] = max(MLlz);
            [~,MLLL] = min(MLlz);
            MLlung(MLULLML).SetRegion('lobe','B');
            MLlung(MLULLML).SetRegion('name','LeftIntermedius');
            MLlung(MLLL).SetRegionDescendants('lobe','LLL');

            % identify upper lobe and lingular
            MLULLML2 = MLlung(MLULLML).children;
            [~, ~, MLl2z] = obj.I2S(ClinicalAirways.SkelEnds(MLULLML2));
            % set LML by lowest z descendant
            [~,MLML] = min(MLl2z);
            MLULLML2(MLML).SetRegionDescendants('lobe','LML');
            
            % set remaining descendants to LUL
            MLULLML2(MLML) = [];
            for idx = 1:length(MLULLML2)
%                 [~,MLUL] = max(MLl2z);
                MLULLML2(idx).SetRegionDescendants('lobe','LUL');
            end

            % identify right upper lobe
            MRlung = MB(rightMBI).children;
            [~, ~, MRz] = obj.I2S(ClinicalAirways.SkelEnds(MRlung));
            [~,MRULi] = max(MRz);
            MRlung(MRULi).SetRegionDescendants('lobe','RUL');


            % check ratio of right major bronchi to the RUL bronchus.
            ratio = MRlung(MRULi).stats.euclength/MB(rightMBI).stats.euclength;
            if ratio < 0.5
                warning(['Relatively short right major bronchi detected, ' ...
                    'this case may have abnormal RUL branching, ' ...
                    'check the lobe classification.' ...
                    'No changes have been made.'])
            end

            % subgraph remaining and get endtubes
            [~,MRMLRLLi] = min(MRz);
            subtubes = MRlung(MRMLRLLi).Descendants;
            endtubes = subtubes(cellfun(@isempty, {subtubes.children}));
            % compare endpoints of endtubes
            endtubes_ep = ClinicalAirways.SkelEnds(endtubes);
            [~, epy, epz] = obj.I2S(endtubes_ep);
            z_minus_y = epz - epy;
            [~, RMLepi] = max(z_minus_y);
            RML_eptube = endtubes(RMLepi);
            [~, RLLepi] = min(z_minus_y);
            RLL_eptube = endtubes(RLLepi);
            % get ancestor list and find intersection
            RML_endpath = [RML_eptube.Ancestors().ID];
            RLL_endpath = [RLL_eptube.Ancestors().ID];
            [intersected_B, intersections] = intersect(RML_endpath, RLL_endpath);

            % assign right bronchus intermedius to branches between carina
            % and RML except the major right bronchus
            for id = intersected_B(2:end)
                obj.tubes(id).SetRegion('lobe','B');
                obj.tubes(id).SetRegion('name','RightIntermedius');
            end

            % assign RML
            RML_id = RML_endpath(min(intersections)-1);
            obj.tubes(RML_id).SetRegionDescendants('lobe','RML');

            % assign remaining labels to RLL
            for ii = 1:length(subtubes)
                if ~isfield(subtubes(ii).region,'lobe')
                    subtubes(ii).SetRegion('lobe','RLL');
                end
            end

            % check for empty lobe fields
            for ii = 1:length(obj.tubes)
                if ~isfield(obj.tubes(ii).region,'lobe')
                    warning(['Lobe classification failed as some tubes were ' ...
                        'not assigned a lobe. Please manually correct.'])
                end
            end

            % set gen by lobe
            obj.RunAllTubes('SetRegionGeneration', 'lobe')
        end
        
        function X = GraphLobarLayout(obj, h, g)
            % Order lobes in graph plot to standard layout.
            %
            % :class:`network.TubeNetwork.Plot` usually orders the lobes randomly.
            % this function will manipulate x positions of lobes so that
            % they appear in the standard layout. This helps make plots
            % more homogenous and therefore comparable to each other.
            %
            % .. warning:
            %   This may cause edges to cross.
            %
            % .. todo:
            %   Integrate a cross over correction mechanism.
            %
            %
            % Args:
            %   h = graphical object of graph plot
            %   g = digraph object
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> [h, g] = AQnet.Plot(colour='lobe');
            %   >>> AQnet.GraphLobarLayout(h, g)
            %

            X = h.XData;

            % identify root nodes of lobes
            lobe_ids = GetTubeValues(obj, 'lobe');

            % find root nodes of lungs
            lobes = {'RUL', 'RML', 'RLL', 'LLL', 'LML', 'LUL'};
            lobe_origins_id = zeros(length(lobes),1);
            % tube index origins in same order as lobes
            for ii = 1:length(lobes)
                lobe = lobes{ii};
                % first occurance of that lobe idx will be the root of it.
                lobe_origins_id(ii) = find(strcmp(lobe_ids,lobe),1);
            end
            
            % check if lungs in the correct order.
            % this is done by comparing the parent nodes of the RUL and
            % LLL.

            % get x pos of the two nodes
            RUL_edge = find(g.Edges.ID == lobe_origins_id(1));
            R_node = g.Edges.EndNodes(RUL_edge,1);
            LLL_edge = find(g.Edges.ID == lobe_origins_id(4));
            L_node = g.Edges.EndNodes(LLL_edge,1);
            % right must be less than left
            if X(R_node) > X(L_node) 
                % swap the right and left lungs at these nodes.
                swapnodes(R_node, L_node)
            end
                        
            % focus on ordering right lobes

            all_lobe_node_xpos = get_all_lobe_node_xpos();
            % if RUL not in 1st pos, swap RUL and right BI
            % find the two decendents of the RUL node
            [~,I] = min(all_lobe_node_xpos);
            if I ~= 1
                % identify right BI
                RUL_node = lobe_origins_id(1);
                R_node_suc = successors(g,R_node);
                RBI_node = R_node_suc(R_node_suc ~= RUL_node);
                swapnodes(RBI_node, RUL_node)
            end
            % swap RML and RLL if in wrong order
            if get_lobe_node_xpos(2) > get_lobe_node_xpos(3)
                % swap LML node and LLL node
                swapnodes(lobe_origins_id(2), lobe_origins_id(3))
            end
            % 

            % focus on ordering left lobes
            all_lobe_node_xpos = get_all_lobe_node_xpos();
            left_lobe_node_xpos = all_lobe_node_xpos(4:6);
            % if LLL not in 1st pos, swap LLL and left BI
            [~,I] = min(left_lobe_node_xpos);
            if I ~= 1
                % swap LLL node and left BI node
                LLL_node = lobe_origins_id(4);
                L_node_suc = successors(g,L_node);
                LBI_node = L_node_suc(L_node_suc ~= LLL_node);
                swapnodes(LBI_node, LLL_node)
            end

            % swap LML and LUL if in wrong order
            if get_lobe_node_xpos(5) > get_lobe_node_xpos(6)
                % swap LML node and LLL node
                swapnodes(lobe_origins_id(5), lobe_origins_id(6))
            end

            % set XData
            h.XData = X;

            function out = get_lobe_node_xpos(idx)
                    % get a lobes x position given index in lobe cell array.
                    lobe_edge = find(g.Edges.ID == lobe_origins_id(idx));
                    lobe_node = g.Edges.EndNodes(lobe_edge,2);
                    out = X(lobe_node);
            end

            function out = get_all_lobe_node_xpos()
                % get all lobes x postion
                out = zeros(length(lobes),1);
                for ii = 1:length(lobes)
                    out(ii) = get_lobe_node_xpos(ii);
                end
            end

            function swapnodes(node_idx_1, node_idx_2)
                % swaps nodes positions in x of two nodes and all 
                % subsequent nodes.

                % identify all nodes on 1
                v1 = dfsearch(g,node_idx_1);

                % identify all nodes on 2
                v2 = dfsearch(g,node_idx_2);

                % identify relative shift
                xshift = X(node_idx_1) - X(node_idx_2);

                % execute relative shift.
                X(v1) = X(v1) - xshift;
                X(v2) = X(v2) + xshift;

            end

        end

        function node_table = ExportGraph(obj, path)
        % export the airway network as a graph. If lobe classifications
        % are present, they will be included.
        %

        parse_filename_extension(path, '.csv');

        % construct node table
        id = [obj.tubes.ID]';
        xyz = cell(length(id),1);
        edge = cell(length(id),1);
        for ii = id'
            % get xyz position from midpoint of tube
            skelpoints = obj.tubes(ii).skelpoints;
            midpoint = floor(length(skelpoints)/2);
            point_pix = obj.I2S(skelpoints(midpoint));
            point_mm = point_pix .* obj.voxdim;
            xyz{ii,1} = num2str(point_mm,'% .2f');
            % convert children nodes to string
            try
                children = num2str([obj.tubes(ii).children.ID]);
            catch
                children = '';
            end
            try
                parent = num2str(obj.tubes(ii).parent.ID);
            catch
                parent = '';
            end
            % save edge
            edge{ii,1} = [parent, ' ', children];
        end
        node_table = table(id, edge, xyz);

        try
            % add lobe classifications if exist
            lobe = AirQuant.list_property({obj.tubes.region},'lobe')';
            node_table = addvars(node_table, lobe);
        catch
        end
            
        % export to csv
        writetable(node_table, path)

        end
    end

    methods(Static)
        function skelends = SkelEnds(airwaylist)
            skelends = cell2mat(cellfun(@(c) [c(end)], ...
                {airwaylist.skelpoints}, 'UniformOutput', false));
        end

    end
end


