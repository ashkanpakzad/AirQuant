% Lead Author: Ashkan Pakzad 2022. ashkanpakzad.github.io.
% See https://github.com/ashkanpakzad/AirQuant for more information.

classdef TubeNetwork < handle
    % TubeNetwork
    %
    % TubeNetwork creates and manages objects inherited from 
    % :class:`tube` that represent anatomical tubes on 3 dimensional
    % images.
    %
    % .. note:: 
    %
    %   TubeNetwork class is intended as a base class for analysing
    %   anatomical tubes. sub classes refined for analysing a particular 
    %   anatomy should be used. e.g. :class:`AirwayNetwork`
    %
    % .. todo:: 
    %   * Make segmentation import more generalised by removing need
    %   for full connectivity.
    %   * Classify seg by tube
    %   * Consider making tubes property private
    %
    % Args:
    %   source: source image, e.g. CT
    %   sourceinfo(struct): header information from source
    %   voxdim: `[float, float, float]` voxel dimensions in mm usually
    %   seg: binary airway segmentation in the same grid space as
    %       :attr:`source` dimensions must match with :attr:`source`.
    %   skel: skeleton based on segementation with no internal loops in the
    %       same grid space as :attr:`source`. Dimensions must match with
    %       :attr:`source` .
    %   lims: reduced volume indices
    %   max_plane_sz: max interpolated slice size
    %   plane_sample_sz: interpolated slice pixel size
    %   spline_sample_sz: mm interval to sample along branch
    %       arclength
    %   skel_points: *private*, list of skeleton points
    %   Dmap: float, array *private*, distance transform of seg
    %
    %
    properties
        source
        sourceinfo
        voxdim
        seg
        skel
        lims
        max_plane_sz
        plane_sample_sz
        spline_sample_sz
        min_tube_sz
        tubes
    end
    properties (SetAccess = private)
        skel_points
        Dmap
    end

    methods
        % init
        function obj = TubeNetwork(source, sourceinfo, seg, skel)
            % Initialise the TubeNetwork class object.
            %
            % The network digraph is constructed using the default method.
            %
            % .. todo:: Expose digraph method to user at initialisation.
            %
            % Args:
            %   source (3darray): CT loaded from nifti using niftiread.
            %   sourceinfo (struct): CT metadata loaded from nifti using
            %       niftiinfo.
            %   seg (3darray): Binary airway segmentation in the
            %       same grid space as CT. Dimensions must match with CT.
            %   skel (3darray): Binary airway centreline in the
            %       same grid space as CT. Dimensions must match with CT
            %
            %

            % check inputs
            seg = logical(seg);
            skel = logical(skel);

            assert(ndims(source) == 3, 'source must be a 3D array.')
            assert(ndims(seg) == 3, 'seg must be a 3D array.')
            assert(ndims(skel) == 3, 'skel must be a 3D array.')

            assert(all(~seg,'all') == false, 'seg is all zero.')
            assert(all(~skel,'all') == false, 'skel is all zero.')
            assert(all(size(source)==size(seg)),['Size of seg ',size(seg), ...
                ' differs from source ', size(source)])
            assert(all(size(source)==size(skel)),['Size of skel ', ...
                size(skel),' differs from source ', size(source)])

            obj.sourceinfo = sourceinfo;

            % process segmentation
            robustseg = ParseSeg(seg);

            % reorient volumes and get properties
            [obj.source, obj.voxdim] = ReorientVolume(source, obj.sourceinfo);
            obj.seg = ReorientVolume(robustseg, obj.sourceinfo);
            obj.skel = ReorientVolume(skel, obj.sourceinfo);

            % identify cropped size by seg
            obj.lims = CropVol(obj.seg);

            % crop all images
            obj.seg = (CropVol(obj.seg, obj.lims));
            obj.source = CropVol(obj.source, obj.lims);
            obj.skel = CropVol(obj.skel, obj.lims);

            % Compute distance transform
            obj = MakeDistanceTransform(obj);

            % Set dynamic resampling parameters and limits
            measure_limit = floor((min(obj.voxdim)/2)*10)/10;
            obj.plane_sample_sz = measure_limit;
            obj.spline_sample_sz = measure_limit;
            obj.min_tube_sz = 3*max(obj.voxdim);

            % Convert skel into digraph
            [digraphout, glink, gnode] = MakeDigraph(obj);

            % .. todo:: Classify segmentation by tubes

            % make tube objects
            obj.tubes = cell(length(glink),1);
            for ii = 1:length(glink)
                obj.tubes{ii,1} = Tube(obj, glink(ii).point);
            end

            % set tube relationships
            
            % classify tubes by generation
%             obj = ComputeAirwayGen(obj);

        end
        
        function obj = MakeDistanceTransform(obj)
            % Compute distance transform of segmentation and save as
            % property.
            obj.Dmap = bwdist(~obj.seg);
        end

        function branch_seg = ClassifySegmentation(obj)
            % .. :todo:: consider making this function more robust!
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

        function obj = ComputeSkelPoints(obj)
            [XP, YP, ZP] = ind2sub(size(obj.seg), find(obj.skel == 1));
            obj.skel_points = [XP, YP, ZP]; % list of skel points
        end

        % GRAPH NETWORK 
        
        function [digraphout, glink, gnode] = MakeDigraph(obj, method)
            if nargin < 2
                method = 'topnode';
            end
            [digraphout, glink, gnode] = Skel2Digraph(obj.skel, method);
        end

        function obj = ComputeAirwayGen(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Indexes each airway generation by shortest path from carina
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop through graph nodes
            G = obj.Gdigraph;
            gens = zeros(length(obj.Glink),1);
            gens(obj.trachea_path) = 0;
            for i = 1:height(G.Nodes)
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

        function SkelAngle(obj)
            % make vector for all airways
            allX = zeros(length(obj.Glink),2);
            allY = zeros(length(obj.Glink),2);
            allZ = zeros(length(obj.Glink),2);
            segvec = zeros(length(obj.Glink),3);
            for ii = 1:length(obj.Glink)
                % identify origin and sink for each link by node
                n1 = obj.Glink(ii).n1;
                n2 = obj.Glink(ii).n2;
                allX(ii,:) = [obj.Gnode(n1).comy, obj.Gnode(n2).comy];
                allY(ii,:) = [obj.Gnode(n1).comx, obj.Gnode(n2).comx];
                allZ(ii,:) = [obj.Gnode(n1).comz, obj.Gnode(n2).comz];
                % convert to vector and store
                segvec(ii,:) = [obj.Gnode(n2).comy - obj.Gnode(n1).comy,...
                    obj.Gnode(n2).comx - obj.Gnode(n1).comx, ...
                    obj.Gnode(n2).comz - obj.Gnode(n1).comz];
            end
            anoms = zeros(length(obj.Glink),1);
            theta = zeros(length(obj.Glink),1);
            for ii = 1:length(obj.Glink)
                if any(ii == obj.trachea_path)
                    continue
                end
                % get angle between airway and (opposite) parent
                p = -segvec(obj.Glink(ii).parent_idx,:);
                c = segvec(ii,:);

                % theta = acos((dot(p,c))/(norm(p)*norm(c)));
                theta(ii) = atan2(norm(cross(p,c)),dot(p,c));
                % identify if anomalous angle (e.g. acute)
                anoms(ii) = (theta(ii) < pi/2);

                if anoms(ii) == 1
                    Color = [1 0 0]; % R
                else
                    Color = [0 1 0]; % G
                end

                plot3(allX(ii,:), allY(ii,:), allZ(ii,:), ...
                    'LineWidth' ,2, 'Color',Color);

                hold on

                % arrow
                q = quiver3(allX(ii,1), allY(ii,1), allZ(ii,1), ...
                    segvec(ii,1), segvec(ii,2), segvec(ii,3));
                q.Color = Color;
                q.AutoScaleFactor = 0.5;
                q.MaxHeadSize = 1.5;
            end

            % convert theta to string array
            theta_str = string(rad2deg(theta));
            for ii = 1:length(theta_str)
                theta_str(ii) = sprintf('%0.0f',theta_str(ii));
                theta_str(ii) = string([char(theta_str(ii)), char(176)]);
            end

            % plot text
            text(allX(:,1)+segvec(:,1)/2+1, ...
                allY(:,1)+segvec(:,2)/2+1, ...
                allZ(:,1)+segvec(:,3)/2+1, ...
                theta_str, 'Color', [0, 0, 0.8])

            %%% nodes
            X_node = [obj.Gnode.comy];
            Y_node = [obj.Gnode.comx];
            Z_node = [obj.Gnode.comz];
            nums_node = string(1:length(obj.Gnode));
            plot3(X_node,Y_node,Z_node, 'k.', 'MarkerSize', 18, 'Color', 'k');
            text(X_node+1,Y_node+1,Z_node+1, nums_node, 'Color', [0.8, 0, 0])

            %axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            view(80,0)
            axis vis3d
            % undo matlab display flip
            ax = gca;
            ax.XDir = 'reverse';
            % Return anomalous branches
            disp(find(anoms == 1))


        end

        function obj = CullAirways(obj, nodes)
            % remove airways beyond a given node from airway tree
            % identify all edges from node outwards
            Glink_exclude = [];

            % get all links after given nodes
            for node = nodes
                [~, E] = bfsearch(obj.Gdigraph, node, 'edgetonew');
                Glink_exclude = [Glink_exclude; E(:)];
            end

            % convert from graph idx to AQ idx
            link_exclude = obj.Gdigraph.Edges.Label(Glink_exclude);

            % remove links
            RemoveLink(obj, link_exclude);
        end

        function obj = ReclassLobeBranch(obj, branch_idx, lobelabel)
            obj.Glink(branch_idx).lobe = lobelabel;
        end
        
        function obj = ReclassLobe(obj, node, lobelabel, dryrun)
            % reclassify airways beyond a given node to a particular
            % lobelabel. This method is useful for when the lobe
            % classifcation algorithm fails. It performs a
            % breadth-first-search (BFS) on the digraph starting at the
            % given node, identifies all airways within it and reclassifys
            % them. It also reclassifies the predecessing edge.
            % Set dryrun = 1 if you want a preview of what edge indices
            % will be changed without effecting any change.

            if nargin < 4
                dryrun = 0;
            end
            % identify all edges from node outwards
            [~, E] = bfsearch(obj.Gdigraph, node, 'edgetonew');
            % add predecessing edge
            E = [E; inedges(obj.Gdigraph, node)];
            % convert graph index to glink index
            GE = [obj.Gdigraph.Edges.Label(E)];

            if dryrun == 1
                printedgeidx = num2cell(GE);
                printlabels = {obj.Glink(GE).lobe}';
                disp([printedgeidx, string(printlabels)])
            else
                % reclass given edges
                newlabels = repmat({lobelabel},1,length(GE));
                [obj.Glink(GE).lobe] = newlabels{:};
                ReclassLobeGen(obj);
            end
        end

        function obj = FilterBranches(obj, percent)
            % Exclude airway branches that are too short/long based on
            % extreme percentage specified. It considers based on their arclengths.
            % This method overwrites property Gdigraph and rewrites Glink and Gnode.
            % Reorganising how airways are connected to eachother.

            % First save/restore Gnode/Glink/Gadj/Gdigraph
            % Ensures that every original airway is considered.
            if isempty(obj.OriginalGraphMap)
                SaveOriginalGraphMap(obj);
            else
                RestoreOriginalGraphMap(obj);
            end

            % Identify airways to exclude
            all_lengths = [obj.Glink.tot_arclength];
            link_exclude_bool = isoutlier(all_lengths, ...
                'percentiles',[percent, 100]);

            % stop crucial airways from conversion
            % e.g. branch airways

            % get list of filtered out airways
            idx_list = 1:length(obj.Glink);
            link_exclude = idx_list(link_exclude_bool);

            RemoveLink(obj, link_exclude)

            % re-run generation classification with new tree
            ReclassLobeGen(obj)
        end

        function RemoveLink(obj, link_exclude)
            % removes branches in link_exclude from AirQuant
            if ~isfield(obj.Glink, 'exclude')
                % add exclusion field to Gnode and Glink
                fielddummy = num2cell(zeros(length(obj.Glink),1));
                [obj.Glink(:).exclude] = fielddummy{:};
                fielddummy = num2cell(zeros(length(obj.Gnode),1));
                [obj.Gnode(:).exclude] = fielddummy{:};
            end

            % Ensure that links are processed from top to bottom.
            [~,E] = bfsearch(obj.Gdigraph,1,'edgetonew');
            % convert from G indices to AQ
            AQ_BFS_E = obj.Gdigraph.Edges.Label(E);
            link_exclude_BFS = intersect(AQ_BFS_E,link_exclude,'stable');


            for ii = link_exclude_BFS'
                % register link to exclude
                obj.Glink(ii).exclude = 1;
                % identify n1 and children
                n1_ii = obj.Glink(ii).n1;
                n2_ii = obj.Glink(ii).n2;
                % register node to exclude
                obj.Gnode(n2_ii).exclude = 1;

                parent = obj.Glink(ii).parent_idx;
                children = obj.Glink(ii).child_idx;

                % remove n2 node from n1 conn
                obj.Gnode(n1_ii).conn(...
                    find(obj.Gnode(n1_ii).conn == n2_ii)) = [];

                % remove child idx from parent of current
                obj.Glink(parent).child_idx(...
                    find(obj.Glink(parent).child_idx == ii)) = [];

                % DIGRAPH - identify edge and n2 graph index
                % remove ii and node 2
                ii_edge = find(obj.Gdigraph.Edges.Label == ii);
                ii_node2 = find(obj.Gdigraph.Nodes.label == n2_ii);
                ii_node1 = find(obj.Gdigraph.Nodes.label == n1_ii);


                for child = children
                    % LINKS
                    % link child branches to n1
                    obj.Glink(child).n1 = n1_ii;
                    n2 = obj.Glink(child).n2;

                    % new parent is parent of current branch
                    obj.Glink(child).parent_idx = parent;

                    % write child idx of new parent
                    obj.Glink(parent).child_idx;

                    % NODES
                    % write new node links
                    obj.Gnode(n1_ii).links(end+1) = child;

                    % connect child's n2 to new n1
                    if all(obj.Gnode(n2).conn ~= n1_ii)
                        obj.Gnode(n2).conn(end+1) = n1_ii;
                    end

                    %%% DIGRAPH edges
                    % identify child G idx, child_node_n2
                    child_edge = find(obj.Gdigraph.Edges.Label == child);
                    child_edge_weight = obj.Gdigraph.Edges.Weight(child_edge);
                    child_node_n2 = find(obj.Gdigraph.Nodes.label == n2);
                    % modify endnodes
                    NewEdgeTable = table([ii_node1, child_node_n2],[child_edge_weight],[child],...
                        'VariableNames',{'EndNodes','Weight','Label'});
                    obj.Gdigraph = rmedge(obj.Gdigraph,child_edge);
                    obj.Gdigraph = addedge(obj.Gdigraph,NewEdgeTable);

                end
                % DIGRAPH - remove after children have been modified
                % remove ii and node 2
                obj.Gdigraph = rmedge(obj.Gdigraph,ii_edge);
                obj.Gdigraph = rmnode(obj.Gdigraph,ii_node2);

            end
        end

        % UTILITIES 

        function out = SuccessReport(obj)
            % produce a table showing success of processing for each airway
            % branch.
            report = struct('airway',num2cell(1:length(obj.Glink)));
            notanyall = @(x) ~any(x,'all');
            % check for each airway arclength and FWHM failures
            for i = 1:length(obj.Glink)
                if any(i == obj.trachea_path)
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

        % HIGH LEVEL - group run lower level methods
        % * spline?
        % * perpinterp
        % * measure

        function obj = AirwayImageAll(obj)
            % Traverse all airway segments except the trachea.
            disp('Start traversing all airway segments')
            total_branches = length(obj.Glink);
            % check to see if any branches already processed.
            incomplete = cellfun(@isempty, obj.TraversedImage);

            for i = 1:length(obj.Glink)
                % skip the trachea or already processed branches
                if any(i == obj.trachea_path) || incomplete(i) == 0
                    disp(['Traversing: ', num2str(i), ' trachea skipped or already complete'])
                    % incase trachea is last branch
                    continue
                end
                obj = CreateAirwayImage(obj, i);
                disp(['Traversing: Completed ', num2str(i), ' of ', num2str(total_branches)])
                % save obj to disk after every 50 branches
                if rem(i,50) == 0
                    save(obj)
                    disp(['SAVE CHECKPOINT at ', num2str(i), ' of ', num2str(total_branches)])
                end
            end
            save(obj)
            disp('Traversing: SAVE CHECKPOINT and Done')
        end
        
        % EXTENT ANALAYSIS
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

        % SEGMENTAL ANALYSIS
        % segmental diameter tapering
        function [intrataper, averagediameter] = ComputeIntraTaperAll(obj, prunelength)
            if nargin < 2
                prunelength = [0 0];
            end
            % loop through branches
            intrataper = NaN(length(obj.arclength), 3);
            averagediameter = NaN(length(obj.arclength), 3);
            for ii = 1:length(obj.arclength)
                disp(ii)
                if any(ii == obj.trachea_path)
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
            prune = (al >= prunelength(1) & al <= al(end) - prunelength(2));
            al = al(prune);
            coeff = NaN(2,3);
            % convert area to diameters
            diameters = real(sqrt(areas/pi)*2);
            for jj = 1:3
                Dvec = diameters(prune,jj);
                try % incase no branch left after pruning/too few points
                    % fit bisquare method
                    coeff(:,jj) = robustfit(al, Dvec,'bisquare');
                    % compute intra-branch tapering as percentage
                    intrataper(jj) = -coeff(2,jj)/coeff(1,jj) * 100;
                catch
                    % leave as NaN
                end
                % compute average
                averagediameter(jj) = trimmean(Dvec, 10);
            end

            if plotflag == 1 && ~any(isnan(intrataper))
                titlevec = ["inner"; "peak"; "outer"];
                for jj = 1:3
                    subplot(3,1,jj)
                    plot(al, diameters(prune,jj), 'k.')
                    hold on
                    plot(al, coeff(2,jj)*al + coeff(1,jj),'r')
                    legend('data', 'bisquare fit', 'Location', 'best')
                    xlabel('arclength (mm)')
                    ylabel('Diameter (mm)')
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
                if any(ii == obj.trachea_path)
                    continue
                end
                for jj = 1:3
                    % identify parent by predecessor node
                    parent = obj.Glink(ii).parent_idx;
                    intertaper(ii,jj) = (averagediameter(parent, jj) - averagediameter(ii,jj))...
                        /(averagediameter(parent, jj)) * 100;
                end
            end

        end

        % segmental volume tapering
        function vol = ComputeIntegratedVol(obj, prunelength, idx)
            % Compute the intertapering value for all airways based on
            % integrated volume along airway.

            if nargin < 2
                prunelength = [0 0];
            end

            % loop through each segment
            vol = NaN(1, 3);

            % get arclength
            arcL = obj.arclength{idx, 1};

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
            prune = (arcL >= prunelength(1) & arcL <= arcL(end) - prunelength(2));
            arcL = arcL(prune);
            % convert area to diameters
            for jj = 1:3
                try % incase no branch left after pruning/too few points
                    % vector of areas
                    Avec = areas(prune,jj);
                    al = arcL;
                    % remove any nan measurements and let integration
                    % 'fill' it in.
                    al(isnan(Avec)) = [];
                    Avec(isnan(Avec)) = [];
                    % integrate arclength against diameter to get volume
                    vol(jj) = trapz(al, Avec);
                catch
                    % leave as NaN
                end
            end
        end

        function allvol = ComputeIntegratedVolAll(obj, prunelength)
            % Compute the intertapering value for all airways based on
            % integrated volume along airway.

            if nargin < 2
                prunelength = [0 0];
            end

            % loop through branches
            allvol = NaN(length(obj.arclength), 3);

            for ii = 1:length(obj.arclength)
                if any(ii == obj.trachea_path)
                    continue
                end
                [allvol(ii,:)] = ComputeIntegratedVol(obj, prunelength, ii);
            end
        end

        function vol_intertaper = ComputeInterIntegratedVol(obj, prunelength)
            % Compute the intertapering value for all airways based on
            % integrated volume along airway.

            if nargin < 2
                prunelength = [0 0];
            end

            % loop through branches
            allvol = NaN(length(obj.arclength), 3);
            vol_intertaper = NaN(length(obj.arclength), 3);

            for ii = 1:length(obj.arclength)
                if any(ii == obj.trachea_path)
                    continue
                end
                [allvol(ii,:)] = ComputeIntegratedVol(obj, prunelength, ii);
            end

            for ii = 1:length(obj.arclength)
                if any(ii == obj.trachea_path)
                    continue
                end
                % identify parent by predecessor node
                parent = find([obj.Glink.n2] == obj.Glink(ii).n1);
                vol_intertaper(ii,:) = (allvol(parent, :) - allvol(ii,:))...
                    ./(allvol(parent, :)) * 100;
            end
        end

        % tortuosity Methods
        function [tortuosity, La, Le] = ComputeTortuosity(obj)
            % La = Arclengths; Le = Euclidean lengths
            La = nan(size(obj.arclength));
            Le = nan(size(La));
            % get difference in euclidean coordinates for each branch.
            for ii = 1:length(La)
                if any(ii == obj.trachea_path)
                    continue
                end
                % get total arc-length
                La(ii) = TotalSplineLength(obj.Splines{ii,1});
                % get euclidean distance
                [~, CT_point_1] = AirQuant.ComputeNormal(obj.Splines{ii,1}, ...
                    obj.Splines{ii,2}(1));
                [~, CT_point_end] = AirQuant.ComputeNormal(obj.Splines{ii,1}, ...
                    obj.Splines{ii,2}(end));
                Le(ii) = norm(CT_point_end - CT_point_1);
            end

            % arclength / euclidean length
            tortuosity = La./Le;
        end

        % generate output
        function SegmentTaperResults = SegmentTaperAll(obj, prunelength)
            % high level function to compute the segmental tapering
            % measurement of all airways.

            % compute taper results by segment
            [intrataper, avg] = ComputeIntraTaperAll(obj, prunelength);
            intertaper = ComputeInterTaper(obj, prunelength);
            vol_intertaper = ComputeInterIntegratedVol(obj, prunelength);
            [tortuosity, arc_length, euc_length] = ComputeTortuosity(obj);
            lobar_intertaper = ComputeLobarInterTaper(obj, prunelength);
            vol = ComputeIntegratedVolAll(obj, prunelength);
            %             parent = [obj.Glink.parent_idx]';

            % organise into column headings
            branch = [1:length(obj.Glink)]';

            inner_intra = intrataper(:, 1);
            peak_intra = intrataper(:, 2);
            outer_intra = intrataper(:, 3);

            inner_avg = avg(:, 1);
            peak_avg = avg(:, 2);
            outer_avg = avg(:, 3);

            inner_inter = intertaper(:, 1);
            peak_inter = intertaper(:, 2);
            outer_inter = intertaper(:, 3);

            inner_lobeinter = lobar_intertaper(:, 1);
            peak_lobeinter = lobar_intertaper(:, 2);
            outer_lobeinter = lobar_intertaper(:, 3);

            inner_volinter = vol_intertaper(:, 1);
            peak_volinter = vol_intertaper(:, 2);
            outer_volinter = vol_intertaper(:, 3);

            inner_vol = vol(:,1);
            outer_vol = vol(:,3);

            if ~isempty(obj.lungvol)
                inner_vol_lung_ratio = inner_vol./obj.lungvol;
                outer_vol_lung_ratio = outer_vol./obj.lungvol;
            else
                inner_vol_lung_ratio = NaN(size(inner_vol));
                outer_vol_lung_ratio = NaN(size(outer_vol));
            end

            thickness_avg = outer_avg - inner_avg;

            % convert to table
            SegmentTaperResults = table(branch, inner_intra, peak_intra, ...
                outer_intra, inner_avg, peak_avg, outer_avg, ...
                inner_inter, peak_inter, outer_inter,...
                inner_volinter, peak_volinter, outer_volinter, ...
                inner_lobeinter, peak_lobeinter, outer_lobeinter, ...
                tortuosity, arc_length, euc_length, inner_vol, outer_vol, ...
                inner_vol_lung_ratio, outer_vol_lung_ratio, thickness_avg);

            % add gen info
            SegmentTaperResults.generation = [obj.Glink.generation]';

            % add lobe info if available
            try % only add lobe information if it exists
                SegmentTaperResults.lobe = [obj.Glink.lobe]';
            catch
            end

            % delete excluded branches
            if isfield(obj.Glink,'exclude')
                SegmentTaperResults(logical([obj.Glink.exclude]),:) = [];
            end

            % Save to AQ object
            obj.Specs.SegmentTaperResults = SegmentTaperResults;
        end

        % TAPERING VISUALISATION 
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
        end

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

        % VISUALISATION
        % Airway Strucutral Tree
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
            h.NodeColor = 'r';
            h.EdgeColor = 'k';
            % check for error edges
            [~,~,erroredge] = DebugGraph(obj);

            % highlight by lobe colour if available
            if isfield(obj.Glink, 'lobe')
                [h, G] = SetGraphLobeColourmap(obj, h, G);
            end

            h.LineWidth = 3;

            if ~isempty(erroredge)
                highlight(h,'Edges',erroredge,'EdgeColor','r');
            end


        end

        function [h, G] = SetGraphLobeColourmap(obj, h, G)
            % set the colours of a network graph by lobe.
            % G is the graph object and h is the plot.
            % i.e. h = plot(G)

            % get lobe info
            lobes = [obj.Glink(:).lobe];
            % convert to graph indices
            edgelobe = lobes(G.Edges.Label);
            % convert labels into numbers
            lobeid = {'B','RUL','RML','RLL','LUL','LML','LLL'};
            cdata = zeros(size(edgelobe));
            for ii = 1:length(cdata)
                [~, ~, cdata(ii)] = intersect(edgelobe(ii),lobeid);
            end

            % set edge colour by index
            G.Edges.EdgeColors = cdata';
            h.EdgeCData = G.Edges.EdgeColors;

            % set colours map and text
            clims = [1 max(cdata(:))];
            colorbarstring = 'Lobe';
            colourshow = clims(1):clims(2);
            colourlabels = lobeid;
            maptype = 'qualitative';
            map = linspecer(max(cdata(:)), maptype);
            colormap(map)
            c = colorbar('Ticks', colourshow, 'TickLabels', colourlabels);
            c.Label.String = colorbarstring;
            caxis(clims)
        end

        function PlotTree(obj, gen, show_seg_txt, show_node_txt)
            % Plot the airway tree with nodes and links in image space. Set
            % gen to the maximum number of generations to show.
            % Original Function by Ashkan Pakzad on 27th July 2019.

            % for max generations set gen to 0
            % hide segment and node labels using show_seg_txt and
            % show_node_txt respectively.

            % does not show excluded airways

            % set defaults
            if nargin < 2
                gen = 0;
            end

            if nargin < 3
                show_seg_txt = 1;
            end

            if nargin < 4
                show_node_txt = 1;
            end

            % if gen is 0, update to max gen
            if gen == 0
                gen = max([obj.Glink(:).generation]);
            end

            % set up reduced link graph and skel
            vis_Glink_logical = [obj.Glink(:).generation] <= gen;
            if isfield(obj.Glink,'exclude')
                vis_Glink_exclude = [obj.Glink(:).exclude] ;
            else
                vis_Glink_exclude = zeros(size(vis_Glink_logical));
            end
            vis_Glink_ind = find(vis_Glink_logical == 1 & vis_Glink_exclude == 0);
            vis_Glink = obj.Glink(vis_Glink_ind);
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

            isosurface(vis_skel)
            alpha(0.7)

            hold on

            % edges
            if show_seg_txt ~= 0
                ind = zeros(length(vis_Glink), 1);
                for i = 1:length(vis_Glink)
                    ind(i) = vis_Glink(i).point(ceil(end/2));
                end
                [Y, X, Z] = ind2sub(size(obj.skel),ind);
                nums_link = string(vis_Glink_ind);
                plot3(X,Y,Z, 'b.', 'MarkerFaceColor', 'none');

                text(X+1,Y+1,Z+1, nums_link, 'Color', [0, 0, 0.8])
            end

            % nodes

            X_node = [vis_Gnode.comy];
            Y_node = [vis_Gnode.comx];
            Z_node = [vis_Gnode.comz];
            nums_node = string(vis_Gnode_ind);
            plot3(X_node,Y_node,Z_node, 'r.', 'MarkerSize', 18, 'Color', 'r');
            if show_node_txt ~= 0
                text(X_node+1,Y_node+1,Z_node+1, nums_node, 'Color', [0.8, 0, 0])
            end

            view(80,0)
            axis vis3d
            % undo matlab display flip
            ax = gca;
            ax.XDir = 'reverse';

        end

        function plot3(obj, gen, show_node_txt)
            % Plot the airway tree in graph form, in 3D. nodes are in
            % in image space. Set gen to the maximum number of
            % generations to show.
            % Original Function by Ashkan Pakzad on 27th July 2019.


            if nargin < 2
                gen = 0;
            end

            if nargin < 3
                show_node_txt = 1;
            end

            if gen == 0
                gen = max([obj.Glink(:).generation]);
            end

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

            %%% e3 edges
            % Set-up lobe colours
            lobeid = {'B','RUL','RML','RLL','LUL','LML','LLL'};
            colours = linspecer(length(lobeid), 'qualitative');
            colormap(colours)
            for i = 1:length(vis_Glink)
                % get lobe colour index
                cidx = strcmp(lobeid, obj.Glink(i).lobe);
                % identify origin and sink for each link by node
                n1 = vis_Glink(i).n1;
                n2 = vis_Glink(i).n2;
                % plot line for each link
                X = [vis_Gnode(n1).comy, vis_Gnode(n2).comy];
                Y = [vis_Gnode(n1).comx, vis_Gnode(n2).comx];
                Z = [vis_Gnode(n1).comz, vis_Gnode(n2).comz];
                plot3(X,Y,Z,'LineWidth',2,'Color', colours(cidx,:))
                hold on
                % get arrow
                U = vis_Gnode(n2).comy - vis_Gnode(n1).comy;
                V = vis_Gnode(n2).comx - vis_Gnode(n1).comx;
                W = vis_Gnode(n2).comz - vis_Gnode(n1).comz;
                q = quiver3(X(1),Y(1),Z(1),U,V,W);
                q.Color = colours(cidx,:);
                q.AutoScaleFactor = 0.5;
                q.MaxHeadSize = 1.5;

            end
            % colorbar
            clims = [1 length(lobeid)];
            c = colorbar('Ticks', clims(1):clims(2), 'TickLabels', lobeid);
            c.Label.String = 'Lobe';
            caxis(clims)


            %%% nodes
            X_node = [vis_Gnode.comy];
            Y_node = [vis_Gnode.comx];
            Z_node = [vis_Gnode.comz];
            nums_node = string(vis_Gnode_ind);
            plot3(X_node,Y_node,Z_node, 'r.', 'MarkerSize', 18, 'Color', 'r');
            if show_node_txt == 1
                text(X_node+1,Y_node+1,Z_node+1, nums_node, 'Color', [0.8, 0, 0])
            end
            %axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            view(80,0)
            axis vis3d
            % undo matlab display flip
            ax = gca;
            ax.XDir = 'reverse';

        end

        % Splines
        function PlotSplineTree(obj)
            % loop through every branch, check spline has already been
            % computed, compute if necessary. Skip trachea. Plot spline.
            % .. :warning: This may not work as well when VoxelSize ~= [1,1,1]
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
            view(80,0)
            axis vis3d
            % undo matlab display flip
            ax = gca;
            ax.XDir = 'reverse';

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
        function PlotSeg(obj)
            % plot segmentation
            patch(isosurface(obj.seg),'EdgeColor', 'none','FaceAlpha',0.1, 'LineStyle', 'none');
            vol3daxes(obj)
        end

        function PlotMap3D(obj, mode)
            % Recommend to use View3D if colour labels appear buggy.

            % mode = 'TaperGradient', 'generation', 'lobes'
            %             axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])

            % generating the color data
            cdata = zeros(size(obj.seg));
            branch_seg = ClassifySegmentation(obj);
            switch mode
                case 'tapergradient'
                    % .. :todo:: rewrite this bit....
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
                    lobeid = {'B','RUL','RML','RLL','LUL','LML','LLL'};
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
            p = patch(isosurface(cdata > 0));

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
                    lobeid = {'B','RUL','RML','RLL','LUL','LML','LLL'};
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
            patch(isosurface(obj.seg),'EdgeColor', 'none','FaceAlpha',0.1);
            hold on
            isosurface(obj.skel,'color','c')
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
            axis square


            try % try block incase FWHMesl has not been executed.
                % plot ray cast results
                hold on
                plot(obj.FWHMesl{link_index, 1}{slide, 1}.x_points + min_centre, obj.FWHMesl{link_index, 1}{slide, 1}.y_points + min_centre,'r.')
                %                 plot(obj.FWHMesl{link_index, 2}{slide, 1}.x_points, obj.FWHMesl{link_index, 2}{slide, 1}.y_points,'c.')
                %                     plot(obj.FWHMesl{link_index, 3}{slide, 1}.x_points + min_centre, obj.FWHMesl{link_index, 3}{slide, 1}.y_points + min_centre,'y.')

                % plot ellipse fitting
                ellipse(obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(4),...
                    obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(1)+min_centre,...
                    obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(2)+min_centre,'m');

                %                 ellipse(obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(4),...
                %                     obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(1),...
                %                     obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(2),'b');

                %                     ellipse(obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(4),...
                %                         obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(1)+min_centre,...
                %                         obj.FWHMesl{link_index, 3}{slide,
                %                         1}.elliptical_info(2)+min_centre,'y');s
                %                     plot(obj.FWHMesl{link_index, 4}{slide, 1}(1)+min_centre,...
                %                         obj.FWHMesl{link_index, 4}{slide, 1}(2)+min_centre, ...
                %                         '.g', 'MarkerSize',20)
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
            awyarray = airwaystack(obj,link_index);
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

        function h = ReformatAirway(obj,link_index,slice_idx)
            % get reformatted airway stack
            awyarray = airwaystack(obj,link_index);
            if nargin < 3
                % set default to middle
                slice_idx = round(size(awyarray,1)/2);
            end
            % generate image
            img = squeeze(awyarray(slice_idx,:,:));
            x = [0 obj.arclength{link_index,1}(end)];
            y = [0 obj.max_plane_sz];
            h = imagesc(x, y, img);
            colormap('gray')
        end

        function awyarray = airwaystack(obj,link_index)
            % generate an airway's interpolated slices into an array
            % stack.
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
        end

        %%% Novel/tapering visualisation
        function [h, G] = GraphPlotDiameter(obj, showlabels, XData, YData)
            if nargin < 2
                showlabels = 1;
            end
            if nargin < 4
                XData = [];
                YData = [];
            end
            % graph plot any variable for each airway as desired. i.e.
            % provide var which is a vector the same length as the number
            % of airways.

            G = obj.Gdigraph;

            if ~exist('obj.Specs.SegmentTaperResults', 'var')
                tapertable = SegmentTaperAll(obj, [0 0]);
            else
                tapertable = obj.Specs.SegmentTaperResults;
            end

            % generate corresponding edgelabels
            if showlabels == 1
                edgelabels = [obj.Glink(G.Edges.Label).generation];
            else
                edgelabels = [];
            end
            edgevar = real(tapertable.inner_avg(G.Edges.Label));

            title('Average Inner lumen Diameter')
            if ~isempty(XData) && ~isempty(XData)
                h = plot(G,'EdgeLabel',edgelabels,'XData',XData,'YData',YData);
            else
                h = plot(G,'EdgeLabel',edgelabels, 'Layout', 'layered');
            end
            h.NodeColor = 'r';
            h.EdgeColor = 'k';

            % set linewidth
            edgevar(isnan(edgevar)) = 0.001;
            h.LineWidth = edgevar;

            % highlight by lobe colour if available
            if isfield(obj.Glink, 'lobe')
                [h, G] = SetGraphLobeColourmap(obj, h, G);
            end
        end

        function h = GraphPlot(obj, thickness_var, label_var)
            if nargin < 3
                label_var = 1:length(obj.Glink);
            end
            % graph plot any variable for each airway as desired. i.e.
            % provide var which is a vector the same length as the number
            % of airways.

            G = obj.Gdigraph;

            % generate corresponding edgelabels
            edgelabels = label_var(G.Edges.Label);
            edgevar = thickness_var(G.Edges.Label);

            title('Average Inner lumen Diameter')
            h = plot(G,'EdgeLabel',edgelabels, 'Layout', 'layered');
            h.NodeColor = 'r';
            h.EdgeColor = 'k';

            % set linewidth
            edgevar(isnan(edgevar)) = 0.001;
            h.LineWidth = edgevar;

            % highlight by lobe colour if available
            if isfield(obj.Glink, 'lobe')
                [h, G] = SetGraphLobeColourmap(obj, h, G);
            end
        end

        function h = LobeAvgPlot(obj, metric)
            % visualise average of each generation across upper and lower
            % lobes for intra/inter/avg values.
            if ~exist('obj.Specs.SegmentTaperResults', 'var')
                tapertable = SegmentTaperAll(obj, [0 0]);
            else
                tapertable = obj.Specs.SegmentTaperResults;
            end

            switch metric % 2 = intra, 5 = avg, 8 = inter
                case 'intra'
                    metricidx = 2;
                    ylab = 'Mean Intratapering (%)';

                case 'avg'
                    metricidx = 5;
                    ylab = 'Mean diameter (mm)';

                case 'inter'
                    metricidx = 8;
                    ylab = 'Mean Intertapering (%)';

                otherwise
                    error('choose either intra, avg or inter.')

            end
            % get max number of generations
            maxgen = max(tapertable.generation);

            % set up upper lobe and lower lobe vars
            ulobegenavg = nan(maxgen,1);
            llobegenavg = nan(maxgen,1);

            % loop through table and average per gen per upper/lower region
            uLobelabels = {'RUL', 'LUL'};
            lLobelabels = {'RLL', 'LLL'};

            for ii = 2:maxgen
                currentgen = tapertable(tapertable.generation == ii,:);
                ulobegenavg(ii) = mean(currentgen.(metricidx)(ismember(currentgen.lobe(:), uLobelabels)),'omitnan');
                llobegenavg(ii) = mean(currentgen.(metricidx)(ismember(currentgen.lobe(:), lLobelabels)),'omitnan');
            end

            % get valid gen indices and flip the upperlobe vals
            ulidx = 1:maxgen;
            llidx = 1:maxgen;
            ulidx =  flip(ulidx(~isnan(ulobegenavg)));
            llidx =  llidx(~isnan(llobegenavg));
            ulobegenavg = flip(ulobegenavg);

            lidx = [ulidx,NaN, llidx];
            y = [ulobegenavg(~isnan(ulobegenavg)); NaN; llobegenavg(~isnan(llobegenavg))];

            h = bar(y, 'k'); % plot dropping nans
            ax = gca;

            ax.XTick = 1:length(lidx);
            ax.XTickLabel = [ulidx,0,llidx];
            xlabel('generation')

            xlab1 = 'Upper Lobes';
            xlab2 = 'Lower Lobes';
            annotation('textbox',[0 0 .1 .2],'String',xlab1,'EdgeColor','none')
            annotation('textbox',[.9 0 .1 .2],'String',xlab2,'EdgeColor','none')

            axis tight
            ylabel(ylab);


        end

        function h = GraphPlotIT(obj)

            % graph plot for plotting intertapering relative to
            % edgethicknes.

            G = obj.Gdigraph;
            % get inner intertaper, set edgelabels and thickness
            intertaper = ComputeInterTaper(obj, [0 0]);
            intertaper = intertaper(G.Edges.Label,1);
            edgevar = abs(intertaper)/5;
            edgelabels = round(intertaper);

            title('Intertaper value')
            h = plot(G,'EdgeLabel',edgelabels, 'Layout', 'layered');
            h.NodeColor = 'r';
            h.EdgeColor = 'k';

            % set linewidth
            edgevar(isnan(edgevar)) = 0.001;
            h.LineWidth = edgevar;

            % highlight negative IT vals
            neg_IT = (intertaper < 0);
            T = G.Edges.EndNodes(neg_IT,:);
            highlight(h,T(:,1),T(:,2),'EdgeColor','r');

        end
        
        % EXPORT 
        function obj = SaveAllAwy(obj, mingen, maxgen, prunelength)
            if nargin < 2
                mingen = 0;
            end

            if nargin < 3 || isnan(maxgen)
                maxgen = max([obj.Glink(:).generation]);
            end

            if nargin < 4
                prunelength = [0 0];
            end
            % make directory
            [fPath, saveid, ~] = fileparts(obj.savename);
            dirname = fullfile(fPath,'airway_patches');
            if ~exist(dirname, 'dir')
                mkdir(dirname)
            end

            % loop through each airway seg
            for ii = 1:size(obj.TraversedImage,1)
                seggen = obj.Glink(ii).generation;
                if  seggen <= mingen || seggen >= maxgen
                    continue
                end

                % choose which slices to save
                al = obj.arclength{ii, 1};
                prune = (al >= prunelength(1) & al <= al(end) - prunelength(2));
                allslices = 1:length(obj.TraversedImage{ii, 1});
                chosenslices = allslices(prune);
                % loop through slices
                for k = chosenslices
                    img = int16(obj.TraversedImage{ii,1}{k,1});

                    % save as int16 TIF
                    imgsavename = fullfile(dirname, [ ...
                        saveid, '_', ...
                        'seg_',num2str(ii), ...
                        '_lobe_', char(obj.Glink(ii).lobe), ...
                        '_gen_', num2str(obj.Glink(ii).generation), ...
                        '_slice_',num2str(k), ...
                        '.tif']);

                    imgdata = img;

                    t = Tiff(imgsavename,'w');
                    tagstruct.Compression = Tiff.Compression.None;
                    tagstruct.ImageLength = size(imgdata,1);
                    tagstruct.ImageWidth = size(imgdata,2);
                    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                    tagstruct.SampleFormat = Tiff.SampleFormat.Int; % int
                    tagstruct.BitsPerSample = 16;
                    tagstruct.SamplesPerPixel = 1;
                    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                    tagstruct.Software = 'AirQuant';
                    setTag(t,tagstruct)
                    write(t,imgdata)
                    close(t);
                    if k == 1
                        disp(imgsavename)
                    end
                end
            end
        end


    end
end