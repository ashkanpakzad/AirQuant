% Lead Author: Ashkan Pakzad 2022. ashkanpakzad.github.io.
% See https://github.com/ashkanpakzad/AirQuant for more information.

classdef TubeNetwork < AirQuant & matlab.mixin.SetGet
    % TubeNetwork for managing analysis of a set of tubes in AirQuant.
    %
    % TubeNetwork creates and manages objects inherited from
    % :class:`tube` that represent anatomical tubes on 3 dimensional
    % images.
    %
    % .. note:
    %   TubeNetwork class is intended as a base class for analysing
    %    anatomical tubes. sub classes refined for analysing a particular
    %    anatomy should be used. e.g. :class:`HumanAirways`
    %
    % .. todo::
    %   * Make segmentation import more generalised by removing need
    %       for full connectivity.
    %   * Consider making tubes property protected.
    %
    % Args:
    %   source(3Darray): source image, e.g. CT
    %   sourceinfo(struct): header information from source
    %   voxdim: `[float, float, float]` voxel dimensions
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
        tubes = [];
        sourceinfo
        skel
        lims
        regioncategories
        spline_sample_sz
        max_plane_sz
        plane_sample_sz
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
        function obj = TubeNetwork(source, sourceinfo, seg, skel, options)
            % Initialise the TubeNetwork class object.
            %
            % The network digraph is constructed using the default method.
            %
            % .. todo::
            %   * Expose digraph method to user at initialisation.
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
            arguments
            source (:,:,:)
            sourceinfo struct
            seg (:,:,:) logical
            skel (:,:,:) logical
            options.fillholes logical = 1
            options.largestCC logical = 0
            end

            assert(ndims(seg) == 3, 'seg must be a 3D array.')
            assert(ndims(skel) == 3, 'skel must be a 3D array.')

            assert(all(~seg,'all') == false, 'seg is all zero.')
            assert(all(~skel,'all') == false, 'skel is all zero.')
            assert(all(size(source)==size(seg)),['Size of seg ',size(seg), ...
                ' differs from source ', size(source)])
            assert(all(size(source)==size(skel)),['Size of skel ', ...
                size(skel),' differs from source ', size(source)])

            obj.sourceinfo = sourceinfo;
            robustseg = ParseSeg(seg, fillholes=options.fillholes, ...
                largestCC=options.largestCC);

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

            % Convert skel into digraph
            [~, glink, ~] = Skel2Digraph(obj);

            % make tube objects
            obj.MakeTubes(glink);

            % set tube relationships
            for ii = 1:length(glink)
                child_idx = find(glink(ii).n2 == [glink(:).n1]);
                for iii = child_idx
                    obj.tubes(ii).SetChildren(obj.tubes(iii));
                end
            end

            obj.RunAllTubes('SetGeneration');

            % classify segmentation to tubes
            obj.ClassifySegmentationTubes();
        end

        function obj = MakeTubes(obj, glink)
            for ii = 1:length(glink)
                obj.tubes = [obj.tubes, Tube(obj, glink(ii).point, ii)];
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
                options.type = 'both'
                options.usesegcrop logical = false
                options.method char = 'linear'
            end

            for ii = progress(1:length(obj.tubes), 'Title', 'Making tube patches')

                if strcmp(options.type,'source') || strcmp(options.type,'both')
                    MakePatchSlices(obj.tubes(ii), obj.source, ...
                        type='source', usesegcrop=options.usesegcrop,...
                        method=options.method);
                end

                if strcmp(options.type,'seg') || strcmp(options.type,'both')
                    MakePatchSlices(obj.tubes(ii), obj.seg, type='seg', ...
                        usesegcrop=options.usesegcrop, method=options.method);
                end

            end

        end

        function obj = ClassifySegmentationTubes(obj)
            % create classified skeleton
            classedskel = zeros(size(obj.skel));
            for ii = 1:length(obj.tubes)
                classedskel(obj.tubes(ii).skelpoints) = obj.tubes(ii).ID;
            end

            % get list of everypoint on segmentation
            [XPQ, YPQ, ZPQ] = obj.I2S(find(obj.seg == 1));
            PQ = [XPQ,YPQ,ZPQ];

            % get list of everypoint on skeleton
            [XP, YP, ZP] = obj.I2S(find(obj.skel == 1));
            P = [XP, YP, ZP];

            % find nearest seg point to it on skeleton
            T = delaunayn(P);
            k = dsearchn(P,T,PQ);

            % find that skeleton point's tube assignment
            classedseg = zeros(size(obj.seg));
            for i = 1:length(PQ)
                skelpoint = P(k(i),:);
                skelval = classedskel(skelpoint(1),skelpoint(2),skelpoint(3));
                classedseg(PQ(i,1),PQ(i,2),PQ(i,3)) = skelval;
            end

            % assign segs to tube
            for ii = 1:length(obj.tubes)
                obj.tubes(ii).segpoints = find(classedseg==obj.tubes(ii).ID);
            end

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
            % .. todo::
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

        function obj = RunAllTubes(obj, tubefunc, varargin)
            % Call a `Tube` method to run on all `obj.tubes`.
            %
            % Useful method for calling a :class:`Tube` method to run on
            % all tubes of a :class:`TubeNetwork` object.
            %
            %
            % Args:
            %   tubefunc (char): CT loaded from nifti using niftiread.
            %   varargin : `OPTIONAL` arguments of method :attr:`tubefunc`.
            %
            %

            assert(isa(tubefunc,"char") || isa(tubefunc,"string"), ...
                'tubefunction must be provided as char/string.')
            for ii = 1:length(obj.tubes)
                obj.tubes(ii).(tubefunc)(varargin{:});
            end
        end

        function I = S2I(obj,I1,I2,I3)
            % short desc
            %
            % long desc
            %
            % .. todo:: add docs
            %
            % Args:
            %   x():
            %
            % Return:
            %   y():
            %

            I = S2I3(size(obj.source),I1,I2,I3);
        end

        function [I1,I2,I3] = I2S(obj,I)
            % short desc
            %
            % long desc
            %
            % .. todo:: add docs
            %
            % Args:
            %   x():
            %
            % Return:
            %   y():
            %

            [I1, I2, I3] = I2S3(size(obj.source),I);
            if nargout == 1
                I1 = [I1; I2; I3]';
            end
        end

        % HIGH LEVEL - group run lower level methods
        % * spline
        % * perpinterp

        function Measure(obj, varargin)
            obj.RunAllTubes('Measure', varargin{:});
        end

        % VISUALISATION - utilites

        function [regionkwarg,regionid] = ParseRegion(obj, regionkwarg)
            % visualisation utility method to interpret the `region`
            % keyword argument
            %
            % looks for a default region that already exists in tubes in
            % order for the user to effortlessy colour figures without the
            % need to be explicit.
            %
            if isempty(regionkwarg) && ~strcmp(regionkwarg,'none') && ...
                ~isempty(fieldnames(obj.tubes(1).region))
                regionfn = fieldnames(obj.tubes(1).region);
                regionkwarg = regionfn{1};
            else
                regionkwarg = [];
            end

            if isfield(obj.regioncategories, regionkwarg)
                regionid = obj.regioncategories.(regionkwarg);
            else
                regionid = [];
            end
        end

        function [cdata, rgb] = ColourIndex(obj, regiontouse, regionid, options)
            % visualisation utility method to set colours by `region`.
            %
            % Args:
            %   maptype(char): default = qua. accepts either 'qua' or 'seq'
            %       to specify the maptype.
            %
            %
            %
            %
            %
            %
            %
            %
            arguments
                obj
                regiontouse
                regionid
                options.maptype char {mustBeMember(options.maptype,{'qua','seq'})} = 'qua'
            end
            % get region info
            regionlist = AirQuant.list_property({obj.tubes.region},regiontouse);

            % convert labels into numbers
            if nargin < 3 || isempty(regionid)
                regionid = unique(AirQuant.list_property({obj.tubes.region},regiontouse));
            end

            cdata = zeros(size(regionlist));
            for ii = 1:length(cdata)
                [~, ~, cdata(ii)] = intersect(regionlist(ii),regionid);
            end

            % set colours map and text
            clims = [1 max(cdata(:))];

            colorbarstring = regiontouse;
            colourshow = clims(1):clims(2);
            colourlabels = regionid;
            rgb = linspecer(max(cdata(:)), options.maptype);
            colormap(rgb)
            c = colorbar('Ticks', colourshow, 'TickLabels', colourlabels);
            c.Label.String = colorbarstring;
            caxis(clims)
        end

        function vol3daxes(obj, ax)
            % utility function for 3D volumetric plotting. Sets the aspect
            % ratio according to voxel size and reverses the x axes for LPS
            % viewing.

            if nargin < 2 % current axes if not specified
                ax = gca;
            end
            axis vis3d
            view(-110, 20)
            % aspect ratio
            ax.DataAspectRatio = 1./obj.voxdim;
            grid on

            f = gcf;
            lh = findobj(f,'Type','light');
            if isempty(lh)
                light;
            end
        end

        % VISUALISATION

        function h = Plot(obj, options)
            arguments
                obj
                options.label = 'ID'
                options.weights = []
                options.shownodes = false
                options.region = ''
            end

            ge = TubesAsEdges(obj);

            switch options.label
                case 'ID' % default
                    edgelabels = ge.Edges.ID;
                case {'generation','gen'}
                    edgelabels = [obj.tubes(ge.Edges.ID).generation];
                otherwise
                    edgelabels = options.label;
            end

            assert(numedges(ge) == length(edgelabels), ['inconsistent ' ...
                'number of edge labels,' num2str(length(edgelabels)), ...
                ' expected ', num2str(numedges(ge))]);

            if options.shownodes == true
                nodelabels = 1:numnodes(ge);
            else
                nodelabels = [''];
            end

            h = plot(ge,nodelabel=nodelabels,edgelabel=edgelabels, ...
                layout='layered');

            h.NodeColor = 'k';
            h.EdgeColor = 'k';

            h.LineWidth = 3;

            % set colour to get chosen region.
            % if region none then output none
            % if region
            [options.region, regionid] = obj.ParseRegion(options.region);

            if ~isempty(options.region)
                cdata = ColourIndex(obj, options.region, regionid);
                edgeregion = cdata(ge.Edges.ID);

                % set edge colour by index
                G.Edges.EdgeColors = edgeregion';
                h.EdgeCData = G.Edges.EdgeColors;
            end

        end

        function Plot3(obj, options)
            % Plot the airway tree in graph form, in 3D. nodes are in
            % in image space. Set gen to the maximum number of
            % generations to show.
            %
            % .. todo: add arrows.
            %
            arguments
                obj
                options.gen = max([obj.tubes.generation])
                options.region = ''
            end

            % set up reduced graph
            vis_Glink_logical = [obj.tubes.generation] <= options.gen;
            vis_tubes = obj.tubes(vis_Glink_logical);

            % set colour to get chosen region.
            [options.region, regionid] = obj.ParseRegion(options.region);

            if ~isempty(options.region)
                [cdata, rgb] = ColourIndex(obj, options.region, regionid);
            else
                rgb = [];
            end

            % plot tube
            for tubeii = vis_tubes
                if ~isempty(rgb)
                    tubeii.Plot3(rgb(cdata(tubeii.ID),:));
                else
                    tubeii.Plot3();
                end
                hold on
            end

            obj.vol3daxes()
            hold off
        end

        function Plot3D(obj, options)
            % Plot volume
            %
            % .. note:
            %   plotting multiple patches (>20) significantly reduces
            %   performance. This is why this function tries to collate every
            %   patch that is designated a different colour and then plots
            %   them.
            %
            %
            %

            arguments
                obj
                options.gen = max([obj.tubes.generation]);
                options.alpha = 0.3
                options.color = 'c'
                options.type {mustBeMember(options.type,{'seg','skel'})} = 'seg'
                options.region = ''
            end


            % set up reduced graph
            vis_Glink_logical = [obj.tubes.generation] <= options.gen;
            vis_tubes = obj.tubes(vis_Glink_logical);

            V = zeros(size(obj.(options.type)));

            % set colour to get chosen region.
            [options.region, regionid] = obj.ParseRegion(options.region);

            if ~isempty(options.region)
                [cdata, rgb] = ColourIndex(obj, options.region, regionid);
            else
                rgb = [];
            end

            % gather tube vol points per color
            for tubeii = vis_tubes
                if ~isempty(rgb)
                    V(tubeii.([options.type,'points'])) = cdata(tubeii.ID);
                else
                    V(tubeii.([options.type,'points'])) = 1;
                end
            end

            % plot each color vol

            if ~isempty(rgb)
                for ii = 1:max(cdata)
                U = zeros(size(V));
                U(V==ii) = 1;
                patch(isosurface(U),...
                'FaceAlpha', options.alpha,...
                'FaceColor', rgb(ii,:),...
                'EdgeColor', 'none');
                hold on
                end
            else
                patch(isosurface(V),...
                'FaceAlpha', options.alpha,...
                'FaceColor', options.color,...
                'EdgeColor', 'none');

            end

            obj.vol3daxes()
            hold off
        end

        function PlotSpline(obj, options)
            % Plot the airway tree in
            %
            %
            %
            arguments
                obj
                options.gen = max([obj.tubes.generation])
                options.region = ''
            end

            % set up reduced graph
            vis_Glink_logical = [obj.tubes.generation] <= options.gen;
            vis_tubes = obj.tubes(vis_Glink_logical);

            % set colour to get chosen region.
            [options.region, regionid] = obj.ParseRegion(options.region);

            if ~isempty(options.region)
                [cdata, rgb] = ColourIndex(obj, options.region, regionid);
            else
                rgb = [];
            end

            % plot tube
            for tubeii = vis_tubes
                if ~isempty(rgb)
                    tubeii.PlotSpline(context=false,color=rgb(cdata(tubeii.ID),:));
                else
                    tubeii.PlotSpline(context=false);
                end
                hold on
            end

            obj.vol3daxes()
            hold off
        end

        function s = OrthoView(obj, type)
            % View volume using MATLAB's inbuilt othogonal viewer.
            %
            % Call MATLAB's orthosliceViewer for a property that is a 3D
            % array of this class. All params are inferred from the object
            % properties.
            %
            % .. todo::
            %   * add documentation to this function
            %   * make scrollable by mousewheel
            %   * refactor to :class:`AirQuant` superclass
            %   * investigate saggital view
            %
            % Args:
            %   type(`char`): `OPTIONAL` default = `source` property name
            %   of 3D array object to view.
            %
            % Return:
            %   - **s** (:attr:`orthosliceViewer`): see ___ for more
            %       details.
            %

            arguments
                obj
                type {mustBeMember(type,{'source','seg'})} = 'source'
            end

            % get volume
            volout = ParseVolOut(obj, type=type);
            volout = permute(volout, [2,1,3]);
            volout = flip(volout, 3);

            % display with orthoview
            s = orthosliceViewer(volout, 'DisplayRangeInteraction','off', ...
                'ScaleFactors',obj.voxdim, 'CrosshairLineWidth', 0.3);

        end

        % Data IO

        % function obj = savetube(obj, tube)
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
            % .. todo:: add documentation to this function
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
            % .. todo:: add documentation to this function
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
