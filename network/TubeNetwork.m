% Lead Author: Ashkan Pakzad 2022. ashkanpakzad.github.io.
% See https://github.com/ashkanpakzad/AirQuant for more information.

classdef TubeNetwork < matlab.mixin.SetGet
    % TubeNetwork
    %
    % TubeNetwork creates and manages objects inherited from
    % :class:`tube` that represent anatomical tubes on 3 dimensional
    % images.
    %
    % .. note:
    %   TubeNetwork class is intended as a base class for analysing
    %    anatomical tubes. sub classes refined for analysing a particular
    %    anatomy should be used. e.g. :class:`AirwayNetwork`
    %
    % .. todo:
    %   * Make segmentation import more generalised by removing need
    %       for full connectivity.
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
        tubes = []
        source
        sourceinfo
        voxdim
        seg
        skel
        lims
        spline_sample_sz
        max_plane_sz
        plane_sample_sz
        min_plane_sz
        plane_scaling_sz = 5;
        savename
    end
    properties (SetAccess = private)
        skel_points
        Dmap
        tubemat
        tubepath
    end

    methods
        % init
        function obj = TubeNetwork(source, sourceinfo, seg, skel)
            % Initialise the TubeNetwork class object.
            %
            % The network digraph is constructed using the default method.
            %
            % .. todo:
            %   * Expose digraph method to user at initialisation.
            %   * Classify segmentation by tubes.
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
            obj.max_plane_sz = 40;
            obj.min_plane_sz = 3*max(obj.voxdim);

            % Convert skel into digraph
            [digraphout, glink, gnode] = Skel2Digraph(obj);

            % make tube objects
            for ii = 1:length(glink)
                obj.tubes = [obj.tubes, Tube(obj, glink(ii).point, ii)] ;
            end

            % set tube relationships
            for ii = 1:length(glink)
                child_idx = find(glink(ii).n2 == [glink(:).n1]);
                for iii = child_idx
                    obj.tubes(ii).SetChildren(obj.tubes(iii));
                end
            end

            obj.RunAllTubes('SetGeneration');

        end

        function obj = RunAllTubes(obj, tubefunc, varargin)
            assert(isa(tubefunc,"char") || isa(tubefunc,"string"), ...
                'tubefunction must be provided as char/string.')
            for ii = 1:length(obj.tubes)
                obj.tubes(ii).(tubefunc)(varargin{:})
            end
        end

        function obj = MakeDistanceTransform(obj)
            % Compute distance transform of segmentation and save as
            % property.
            obj.Dmap = bwdist(~obj.seg);
        end

        function obj = MakeTubePatches(obj, options)
            arguments
                obj
                options.makesegpatches logical = false
                options.usesegcrop logical = false
            end
            for ii = 1:length(obj.tubes)
                MakePatchSlices(obj.tubes(ii,1), obj.source, type='source', usesegcrop=options.usesegcrop);
                if options.makesegpatches == true
                    MakePatchSlices(obj.tubes{ii,1}, obj.seg, type='seg', usesegcrop=options.usesegcrop);
                end
            end
        end

        function obj = ComputeSkelPoints(obj)
            [XP, YP, ZP] = ind2sub(size(obj.seg), find(obj.skel == 1));
            obj.skel_points = [XP, YP, ZP]; % list of skel points
        end

        % GRAPH NETWORK
        function [g, glink, gnode] = Skel2Digraph(obj, method)
            if nargin < 2
                method = 'topnode';
            end
            [g, glink, gnode] = skel_2_digraph(obj.skel, method);
        end

        function g = TubesAsEdges(obj)
            % makes the digraph based on tube relationships

            % init edgetable

            gn = TubesAsNodes(obj);

            % add incomming edge to nodes have no incoming edges
            nin = find(indegree(gn) == 0);

            asedges = gn.Edges;
            for nini = nin
                asedges.EndNodes(height(asedges)+1,:) = [max(asedges.EndNodes(:))+1 nini];
            end
            asedges.ID = asedges.EndNodes(:,2);

            g = digraph(asedges);
        end

        function g = TubesAsNodes(obj)
            % make the digraph with tubes as nodes
            %
            %
            % ..todo:
            %   * Scope for improving efficiency, variables that change
            %       size.
            %

            % get all tubes with generation 0
            zerogentubes = find([obj.tubes(:).generation] == 0);

            % set up initial nodes with zero gen tubes
            s = [];
            t = [];

            % DFS from each zero gen tube, creating edgetable.
            for ii = 1:length(zerogentubes)
                edgestosearch = obj.tubes(zerogentubes(ii));
                while ~isempty(edgestosearch)
                    current_tube = edgestosearch(1);
                    s_parent = current_tube.ID;
                    if ~isempty(current_tube.children)
                        for childtube = current_tube.children
                            s = [s, s_parent];
                            t = [t, childtube.ID];
                            edgestosearch = [edgestosearch, childtube];
                        end
                    end
                    edgestosearch(1) = [];
                end
            end

            g = digraph(s,t);

        end

        % UTILITIES
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

        % HIGH LEVEL - group run lower level methods
        % * spline
        % * perpinterp
        % * measure

        % VISUALISATION
        % Airway Strucutral Tree
        function h = plot(obj, options)
            arguments
            obj
            options.label = 'ID'
            options.weights = []
            options.shownodes = false
            end
            % Default plot is a graph network representation.

            ge = TubesAsEdges(obj);

            edgelabels = ge.Edges.ID;

            switch options.label
                case isnumeric(options.label)
                    edgelabels = options.label;
                case 'ID' % default

                case {'generation','gen'}
                    edgelabels = [obj.tubes(edgelabels).generation];
                otherwise
                    warning('Unexpected options.label type. Resorting to default type.')
            end

            if options.shownodes == true
                nodelabels = 1:numnodes(ge);
            else
                nodelabels = [''];
            end

            h = plot(ge,nodelabel=nodelabels,edgelabel=edgelabels,layout='layered');

            h.NodeColor = 'r';
            h.EdgeColor = 'k';
            h.LineWidth = 3;
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

        % Data IO

%         function obj = savetube(obj, tube)
%             % update specific tube object only in matfile
%             %
%             % firstrun
%             if isempty(obj.tubepath)
%               % init matfile
%               obj.tubemat = matfile(fileparts(obj.objpath,'tubes.mat'),'Writable',true);
%               % save each tube to new variable in mat by ID
%               for ii = 1:length(obj.tubes)
%                 currenttube = obj.tubes(ii);
%                 currentID = ['tube_',num2str(currenttube.ID)];
%                 assignin('base',currentID, currenttube)
%                 obj.tubemat
%               end
%             end
% 
%         end

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

        function save(obj, path)
            % save class object to disk
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %   * consider copying object and removing tubes property.
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            tubes = obj.tubes;
            obj.tubes = [];
            save(path, 'obj', '-v7.3')
            disp(['Saved to ',path])
            % tubes save seperately
            obj.tubepath = replace(path,'.m','_tubes.m');
            save(obj.tubepath, 'tubes');
            obj.tubemat = matfile(obj.tubepath,'Writable',true);
            disp(['Tubes cache saved to ',obj.tubepath]);
            obj.tubes = tubes;
        end
          
        function savetube(obj,tubetosave)
            % change only tube on disk.
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %   * Needs attention.: make tubes individual files to save. 
            %
            %
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %

            assert(~isempty(obj.tubemat), 'Need to first run save methods, see documentation.')
            objid = find([obj.tubes(:).ID] == tubetosave.ID);
            assert(~isempty(objid), 'tube ID not found in network object.')
            obj.tubemat.tubes(1,objid) = tubetosave;
            
        end
    end

    methods(Static)
        function obj = load(path)
          s = load(path);
          obj = s.obj;
          disp(['Network loaded from ', path]);

          obj.tubepath = replace(path,'.m','_tubes.m');
          t = load(obj.tubepath);
          obj.tubemat = matfile(obj.tubepath,'Writable',true);
          obj.tubes = t.tubes;
          disp(['Tubes loaded from ',obj.tubepath]);
        end
    end
end
