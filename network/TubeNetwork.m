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
    %    anatomy should be used. e.g. :class:`ClinicalAirways`
    %
    % .. todo::
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
            % Optional arguments for initialising relating to segmentation
            % processing and then analysis for methods called later. The
            % :attr:`spline_sample_sz` and :attr:`plane_sample_sz`
            % variables are by default set to half the shortest dimension
            % of the source voxel size.
            %
            %
            %
            % Args:
            %   source (3darray): CT loaded from nifti using niftiread.
            %   sourceinfo (struct): CT metadata loaded from nifti using
            %       niftiinfo.
            %   seg (3darray): Binary airway segmentation in the
            %       same grid space as CT. Dimensions must match with CT.
            %   skel (3darray): Binary airway centreline in the
            %       same grid space as CT. Dimensions must match with CT
            %   fillholes (bool): *OPTIONAL* `default = true` whether to
            %       fill any holes found in a 2.5D (any 2D orthogonal plane) 
            %       search.
            %   largestCC (bool): *OPTIONAL* `default = false` Only keep
            %       the largest connected component in the segmentation.
            %   spline_sample_sz (float): *OPTIONAL* `default =
            %   floor((min(obj.voxdim))*10)/10` spline interval to
            %       use for all tubes.
            %   plane_sample_sz (float): *OPTIONAL* `default =
            %   floor((min(obj.voxdim))*10)/10` patch interpolation size
            %       to use for all tubes.
            %   max_plane_sz (float): *OPTIONAL* `default =40` max patch 
            %       size when tubes interpolated in same units as voxel
            %       dimensions.
            %   originmethod (float): *OPTIONAL* `default =
            %       'topnode'` method to set the tube network's origin. Which
            %       will be treated as generation 0. See
            %       :class:`network.TubeNetwork.Skel2Digraph` for options.
            %
            %
            arguments
            source (:,:,:)
            sourceinfo struct
            seg (:,:,:) logical
            skel (:,:,:) logical
            options.fillholes logical = 1
            options.largestCC logical = 0
            options.spline_sample_sz = nan
            options.plane_sample_sz = nan
            options.max_plane_sz = nan
            options.reorient logical = 1
            options.voxdim = nan
            options.originmethod = 'topnode';
            end

            assert(ndims(seg) == 3, 'seg must be a 3D array.')
            assert(ndims(skel) == 3, 'skel must be a 3D array.')

            assert(all(~seg,'all') == false, 'seg is all zero.')
            assert(all(~skel,'all') == false, 'skel is all zero.')
            assert(all(size(source)==size(seg)),['Size of seg ',num2str(size(seg)), ...
                ' differs from source ', num2str(size(source))])
            assert(all(size(source)==size(skel)),['Size of skel ', ...
                size(skel),' differs from source ', num2str(size(source))])

            obj.sourceinfo = sourceinfo;
            
            % process segmentation
            robustseg = ParseSeg(seg, fillholes=options.fillholes, ...
                largestCC=options.largestCC);
            
            % process skeleton
            robustskel = ParseSeg(skel, fillholes=false, ...
                largestCC=options.largestCC);

            % reorient volumes and get properties
            if options.reorient == true
                [obj.source, obj.voxdim] = ReorientVolume(source, obj.sourceinfo);
                obj.seg = ReorientVolume(robustseg, obj.sourceinfo);
                obj.skel = ReorientVolume(robustskel, obj.sourceinfo);
            else
                obj.source = source;
                obj.seg = robustseg;
                obj.skel = robustskel;
                obj.voxdim = sourceinfo.PixelDimensions;
            end
            
            % manual voxdim
            if ~isnan(options.voxdim)
                assert(all(size(options.voxdim)==[1,3]), ...
                    strcat('Size of voxdim must be (1,3), got', ...
                    num2str(size(options.voxdim))))
                obj.voxdim = options.voxdim;
            end
            
            disp('Dimensions in physical units, usually millimetres.')
            disp(['Voxel dimensions = ', ...
                num2str(obj.voxdim)])

            % identify cropped size by seg
            obj.lims = CropVol(obj.seg);

            % crop all images
            obj.seg = (CropVol(obj.seg, obj.lims));
            obj.source = CropVol(obj.source, obj.lims);
            obj.skel = CropVol(obj.skel, obj.lims);

            % Compute distance transform
            obj = MakeDistanceTransform(obj);

            % Set dynamic resampling parameters and limits
            measure_limit = floor((min(obj.voxdim))*10)/10;

            % patch sample size
            if ~isnan(options.plane_sample_sz)
                obj.plane_sample_sz = options.plane_sample_sz;
            else
                obj.plane_sample_sz = measure_limit;
            end
            disp(['Patch sampling size = ', ...
                num2str(obj.plane_sample_sz)])

            % spline sample size
            if ~isnan(options.spline_sample_sz)
                obj.spline_sample_sz = options.spline_sample_sz;
            else
                obj.spline_sample_sz = measure_limit;
            end
            disp(['Spline sampling size = ', ...
                num2str(obj.spline_sample_sz)])
        
            if min(obj.voxdim) < obj.plane_sample_sz || min(obj.voxdim) < obj.spline_sample_sz
                warning(['The smallest voxel diameter is ', ...
                num2str(min(obj.voxdim)), '. This is smaller than the ' ...
                    'patch/spline sample size. This could have unexpected results.'])
            end
            % max patch sample size
            if ~isnan(options.max_plane_sz)
                obj.max_plane_sz = options.max_plane_sz;
            else
                obj.max_plane_sz = 40;
            end
            disp(['Max plane size = ', ...
                num2str(obj.max_plane_sz)])

            % Convert skel into digraph
            [~, glink, ~] = Skel2Digraph(obj, options.originmethod);

            % make tube objects
            obj.MakeTubes(glink);
            
            disp(strcat(num2str(length(glink)), ' tubes found.'))

            % set tube relationships
            for ii = 1:length(glink)
                child_idx = find(glink(ii).n2 == [glink(:).n1]);
                for iii = child_idx
                    obj.tubes(ii).SetChildren(obj.tubes(iii));
                end
            end

            obj.RunAllTubes('SetGeneration');

            % Compute angles of tubes
            obj.RunAllTubes('ComputeDirections');
            obj.RunAllTubes('ComputeChangeAngle');
            obj.RunAllTubes('ComputeParentAngle');
            obj.RunAllTubes('ComputeSiblingAngle',0);

            % classify segmentation to tubes
            obj.ClassifySegmentationTubes();
        end

        function obj = MakeTubes(obj, glink)
            % Internal method for making tubes.
            %
            % This method should be changed for new anatomies to call
            % different types of children of :class:`Tube.tube`.
            %
            for ii = 1:length(glink)
                obj.tubes = [obj.tubes, Tube(obj, glink(ii).point, ii)];
            end
        end

        function obj = MakeDistanceTransform(obj)
            % Compute distance transform of segmentation and save as
            % property :attr:`Dmap`.
            %
            obj.Dmap = bwdist(~obj.seg);
        end

        function obj = MakeTubePatches(obj, options)
            % Calls :meth:`tube.Tube.MakePatchSlices` on all objects in the 
            % :attr:`tubes`.
            %
            % Args:
            %   type(char): *OPTIONAL* `default = 'both'` either `'both'`, 
            %       `'source'` or `'seg'`.
            %   usesegcrop(logical): *OPTIONAL* `default = false` use 
            %       dynamic slicing, vary the size of the plane dependant 
            %       on :attr:`tube.Tube.network.seg` size at that point.
            %   method(char): *OPTIONAL* `default = 'cubic'` see `interp3`_ for details.
            %
            % .. _interp3: https://www.mathworks.com/help/matlab/ref/interp3.html
            %
            arguments
                obj
                options.type = 'both'
                options.usesegcrop logical = false
                options.method char = 'linear'
                options.gpu logical = 1
            end

            for ii = progress(1:length(obj.tubes), 'Title', 'Making tube patches')

                if strcmp(options.type,'source') || strcmp(options.type,'both')
                    MakePatchSlices(obj.tubes(ii), obj.source, ...
                        type='source', usesegcrop=options.usesegcrop,...
                        method=options.method, gpu=options.gpu);
                end

                if strcmp(options.type,'seg') || strcmp(options.type,'both')
                    MakePatchSlices(obj.tubes(ii), obj.seg, type='seg', ...
                        usesegcrop=options.usesegcrop, method=options.method, ...
                        gpu=options.gpu);
                end

            end

        end

        function obj = ClassifySegmentationTubes(obj)
            % Classify :attr:`seg` into its tube components.
            %
            % Classifies every foreground point of :attr:`seg` to the
            % nearest point in :attr:`skel`. Assigning it to the respective
            % :class:`tube.Tube` object in :attr:`tubes`.
            %

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

        function ge = TubesAsEdges(obj)
            % makes the digraph based on tube relationships

            % init edgetable
            gn = TubesAsNodes(obj);
            asedges = gn.Edges;

            % add incomming edge to nodes that have no incoming edges
            nin = find(indegree(gn) == 0);
            
            for ii = 1:length(nin)
                % make new node
                gn = addnode(gn,1);
                % make edge into node
                gn = addedge(gn,numnodes(gn),nin(ii));
%                 asedges.EndNodes(height(asedges)+1,:) = [max(asedges.EndNodes(:))+1 nini];
            end
            
            gn.Edges.ID = gn.Edges.EndNodes(:,2);

            ge = gn;

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
            g = digraph();
            
            % DFS from each zero gen tube, creating edgetable.
            for ii = 1:length(zerogentubes)
                edgestosearch = obj.tubes(zerogentubes(ii));
                if isempty(edgestosearch.parent) && isempty(edgestosearch.children)
                    g = addnode(g,1);
                end
                while ~isempty(edgestosearch)
                    current_tube = edgestosearch(1);
                    parent = current_tube.ID;
                    if ~isempty(current_tube.children)
                        for childtube = current_tube.children
                            g = addedge(g, parent, childtube.ID);
                            edgestosearch = [edgestosearch, childtube];
                        end
                    end
                    edgestosearch(1) = [];
                end
            end
        end

        % UTILITIES

        function obj = RunAllTubes(obj, tubefunc, varargin)
            % Call a :class:`tube.Tube` method to run on all :attr:`tubes`.
            %
            % Provides progress of for-loop running the method.
            %
            % Args:
            %   tubefunc (char): Any method name of :attr:`tubes` object.
            %   varargin: arguments of called method in `tubefunc`.
            %
            % .. todo: add example.
            %

            assert(isa(tubefunc,"char") || isa(tubefunc,"string"), ...
                'tubefunction must be provided as char/string.')
            for ii = progress(1:length(obj.tubes), 'Title', strcat('RunAllTubes: ', tubefunc))
                obj.tubes(ii).(tubefunc)(varargin{:});
            end
        end

        function I = S2I(obj,I1,I2,I3)
            % Fast specific implementation of MATLAB's `sub2ind
            % <https://uk.mathworks.com/help/matlab/ref/sub2ind.html>`_.
            %

            I = S2I3(size(obj.source),I1,I2,I3);
        end

        function [I1,I2,I3] = I2S(obj,I)
            % Fast specific implementation of MATLAB's `ind2sub
            % <https://uk.mathworks.com/help/matlab/ref/ind2sub.html>`_.
            %

            [I1, I2, I3] = I2S3(size(obj.source),I);
            if nargout == 1
                I1 = [I1; I2; I3]';
            end
        end

        function Measure(obj, varargin)
            % run chosen measure class on all :attr:`tubes`. See method
            % :meth:`tubes.Tube.Measure`.
            %
            % Also calls measurement methods that require all tube
            % measurement methods to be completed to run.
            %
            obj.RunAllTubes('Measure', varargin{:});
            obj.RunAllTubes('ComputeIntertaper');
        end

        % VISUALISATION - utilites
        
        function outlabel = GetTubeValues(obj, labelname, labelidx)
            % get the value of each tube property or stat
            
            % search for label in tube properties
            if isprop(obj.tubes(1), labelname)
                object_eval = ('obj.tubes(ii)');
                % then search for label in tube stats
            elseif isfield(obj.tubes(1).stats, labelname)
                object_eval = ('obj.tubes(ii).stats');
                % then search for label in region
            elseif isfield(obj.tubes(1).region, labelname)
                object_eval = 'obj.tubes(ii).region';
            elseif isempty(labelname)
                % if blank then return blank.
                outlabel = {};
                return
            else
                error(['No tube property/stats/region found with name ', labelname,'.'])
            end
            
            % get num tubes and reserve memory
            numtubes = length(obj.tubes);
            outlabel = cell(size(obj.tubes));

            % get values
            for ii = 1:numtubes
                label_all = [eval(object_eval).(labelname)];
                if ischar(label_all)
                    outlabel{ii} = label_all;
                else
                    outlabel{ii} = label_all(labelidx);
                end
            end
            
            % if value type is numerical then convert to matrix
            if isnumeric(outlabel{1})
                outlabel = cell2mat(outlabel);
            end

        end

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

        function [cdata, rgb] = ColourIndex(obj, vals, options)
            % Set colours properties per tube segment.
            %
            % Visualisation utility method to set colours properties per
            % tube segment. It behaves differently depending on
            % quantitative or qualitative input. i.e. change in colourmap
            % style and colourbar tick value style.
            %
            % ..todo:
            %   * Reinstate set order for colours. e.g. region properties.
            %
            % Args:
            %   vals(array or cell array): list of values per tube.
            %   types(array or cell array): *OPTIONAL* `default = ''` 
            %       possible values that `vals` can take. Leave empty to
            %       infer.
            %   name(char): *OPTIOANL* `default = ''` name to appear on
            %       colourbar label.
            %   maptype(char): *OPTIONAL* `default = ''` accepts either
            %       `''`, `'qua'` or `'seq'`. Leave blank to infer by `vals` 
            %       input.
            %
            arguments
                obj
                vals (1,:)
                options.types = ''
                options.name = ''
                options.maptype char {mustBeMember(options.maptype,{'qua','seq',''})} = ''
            end

            assert(length(obj.tubes) == length(vals), ['Should have same number ' ...
                'of vals input as there are tubes. Got ', num2str(length(vals)), ...
                ', but expected ', num2str(length(obj.tubes))])
            
            % convert labels into numbers
            if isempty(options.types)
                types = unique(vals);
            else
                types = options.types;
            end
            
            % bin into unique bins
            cdata = zeros(size(vals));
            for ii = 1:length(cdata)
                [~, ~, cdata(ii)] = intersect(vals(ii),types);
            end

            % infer colourmap type if necessary
            if isempty(options.maptype)
                if iscell(vals)
                    options.maptype = 'qua';
                else
                    options.maptype = 'seq';
                end
            end

            % set colours map and text
            clims = [1 max(cdata(:))];
            colorbarstring = options.name;
            rgb = linspecer(max(cdata(:)), options.maptype);
            colormap(rgb)
            if strcmp(options.maptype, 'qua')
                colourshow = clims(1):clims(2);
                colourlabels = types;
            else
                colourshow = [clims(1), clims(2)];
                colourlabels = [min(types), max(types)];
            end
            c = colorbar('Ticks', colourshow, 'TickLabels', colourlabels);
            c.Label.String = colorbarstring;
            caxis(clims)
        end

        function vol3daxes(obj, ax)
            % utility function for 3D volumetric plotting. Sets the aspect
            % ratio according to voxel size and reverses the x axes for LPS
            % viewing.
            %
            % .. todo: documentation is a stub
            %
            if nargin < 2 % current axes if not specified
                ax = gca;
            end
            axis vis3d
            view(80, 20)
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
            % Very powerful graph visualisation tool of tubes. Can
            % manipulate by label, edge weight and colour.
            %
            % Args:
            %   shownodes(bool) = *OPTIONAL* `default = false`
            %   label = *OPTIONAL* `default = 'ID'` set edge labels.
            %       if `char` must be an :class:`tube` property, :class:`tube` `stats` or `region` field.
            %       e.g. `'generation'`, `'arclength'`, `'lobe'`.
            %       Can also be vector of length equal to number of tubes
            %       in order. 
            %   labelidx(scalar) = *OPTIONAL* `default = 1`. Index of
            %       chosen property in `label`.
            %   weight(char) = *OPTIONAL* `default = none`. set line 
            %       thickness.
            %       if `char` must be an :class:`tube` property, :class:`tube` `stats` or `region` field.
            %       e.g. `'generation'`, `'arclength'`, `'lobe'`.
            %       Can also be vector of length equal to number of tubes
            %       in order. 
            %   weightidx(scalar) = *OPTIONAL* `default = 1`. Index of
            %       chosen property in `weight`.
            %   weightfactor(float) = *OPTIONAL* `default = 1`
            %     determines the highest scaling of the linethickness.
            %   colour(char) = *OPTIONAL* `default = 1`. Set variable for
            %       colour labelling.
            %       if `char` must be an :class:`tube` property, :class:`tube` `stats` or `region` field.
            %       e.g. `'generation'`, `'arclength'`, `'lobe'`.
            %       Can also be vector of length equal to number of tubes
            %       in order. 
            %   colouridx(scalar) = *OPTIONAL* `default = 1` Index of
            %       chosen property in `colour`.          
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.Plot();
            %
            % .. |network_Network_Plot| image:: figs/network_plot.png
            %    :width: 400
            %    :alt: figure plot - Network Plot
            %
            % |network_Network_Plot|
            %
            arguments
                obj
                options.shownodes = false
                options.label = 'ID'
                options.labelidx = 1
                options.weight = ''
                options.weightidx = 1
                options.weightfactor = 1
                options.colour = ''
                options.colouridx = 1
            end
            
            % get graph layout
            ge = TubesAsEdges(obj);

            % get label list
            outlabel = GetTubeValues(obj, options.label, options.labelidx);
            if isempty(outlabel)
                edgelabels = outlabel;
            else
                edgelabels = outlabel(ge.Edges.ID);
                assert(numedges(ge) == length(edgelabels), ['inconsistent ' ...
                    'number of edge labels,' num2str(length(edgelabels)), ...
                    ' expected ', num2str(numedges(ge))]);
            end

            % node labels
            if options.shownodes == true
                nodelabels = 1:numnodes(ge);
            else
                nodelabels = [''];
            end

            % plot
            h = plot(ge,nodelabel=nodelabels,edgelabel=edgelabels, ...
                layout='layered');

            % node colours
            h.NodeColor = 'k';
            h.EdgeColor = 'k';

            outweight = GetTubeValues(obj, options.weight, options.weightidx);
            if isempty(outweight)
                edgevar = 1;
            else
                edgevar = outweight(ge.Edges.ID);
                assert(numedges(ge) == length(edgevar), ['inconsistent ' ...
                    'number of edge weights,' num2str(length(edgevar)), ...
                    ' expected ', num2str(numedges(ge))]);
            end

            % scale up thickest line
            max_thick = max(edgevar);
            scale = options.weightfactor/max_thick;
            edgevar = edgevar*scale;

            % set nan or 0 variable edges to very small value
            edgevar(isnan(edgevar)) = 0.001;
            edgevar(edgevar==0) = 0.001;

            h.LineWidth = edgevar;

            % set colour

            if ~isempty(options.colour)
                try
                    colours = GetTubeValues(obj, options.colour, options.colouridx);
                    colourvals = colours(ge.Edges.ID);
                    cdata = ColourIndex(obj, colourvals, name=options.colour);

                    % set edge colour by index
                    G.Edges.EdgeColors = cdata';
                    h.EdgeCData = G.Edges.EdgeColors;
                catch
                    warning(['Attempted to colour image but failed. Likely ' ...
                        'due to incomplete definition for all tubes.'])
                end
            end

        end

        function Plot3(obj, options)
            % Plot the airway tree in graph form, in 3D. nodes are in
            % in image space. Set gen to the maximum number of
            % generations to show.
            %
            % .. todo: 
            %   * add arrows.
            %   * add weight parameter.
            %   * documentation stub
            %
            % Args:
            %   gen(int): *OPTIONAL* `default =
            %       max([obj.tubes.generation])` plot up to which generation.
            %       All by default.
            %   colour(char) = *OPTIONAL* `default = 1`. Set variable for
            %       colour labelling.
            %       if `char` must be an :class:`tube` property, :class:`tube` `stats` or `region` field.
            %       e.g. `'generation'`, `'arclength'`, `'lobe'`.
            %       Can also be vector of length equal to number of tubes
            %       in order. 
            %   colouridx(scalar) = *OPTIONAL* `default = 1` Index of
            %       chosen property in `colour`.      
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.Plot3();
            %
            % .. |network_Network_Plot3| image:: figs/network_plot3.png
            %    :width: 400
            %    :alt: figure plot - Network Plot3
            %
            % |network_Network_Plot3|
            %
            arguments
                obj
                options.gen = max([obj.tubes.generation])
                options.colour = ''
                options.colouridx = 1
            end

            % set up reduced graph
            vis_Glink_logical = [obj.tubes.generation] <= options.gen;
            vis_tubes = obj.tubes(vis_Glink_logical);

            % set colour

            if ~isempty(options.colour)
                try
                colourvals = GetTubeValues(obj, options.colour, options.colouridx);
                [cdata, rgb] = ColourIndex(obj, colourvals, name=options.colour);
                hold on
                catch
                    warning(['Attempted to colour image. Failed ' ...
                        'due to incomplete region definition for all tubes.'])
                    rgb = [];
                end
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
            % Plot segmentation surface of all :attr:`tubes`.
            %
            %
            % Args:
            %  gen(int): *OPTIONAL* `default =
            %       max([obj.tubes.generation])` plot up to which generation.
            %       All by default.
            %  alpha (float): *OPTIONAL* `default = 0.3` opacity of surface
            %   plot.
            %  colour: *OPTIONAL* `default = 'c'` color of surface. Can be any
            %   MATLAB accepted color format, e.g. RGB.
            %  type(char): *OPTIONAL* `default = 'seg'` can be either
            %  `'seg'` or `'skel'`. which to plot surface.
            %  region
            %
            % .. note:
            %   plotting multiple patches (>20) significantly reduces
            %   performance. This is why this function tries to collate everys
            %   patch that is designated a different colour and then plots
            %   them. 
            % 
            % .. warning:
            %   Choosing a colour option that is continuous will cause this
            %   method to take a very long time.
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.Plot3D();
            %
            % .. |network_Network_Plot3D| image:: figs/network_plot3d.png
            %    :width: 400
            %    :alt: figure plot - Network Plot3D
            %
            % |network_Network_Plot3D|
            %

            arguments
                obj
                options.gen = max([obj.tubes.generation]);
                options.alpha = 0.3
                options.colour = ''
                options.colouridx = 1
                options.type {mustBeMember(options.type,{'seg','skel'})} = 'seg'
            end


            % set up reduced graph
            vis_Glink_logical = [obj.tubes.generation] <= options.gen;
            vis_tubes = obj.tubes(vis_Glink_logical);

            V = zeros(size(obj.(options.type)));

            % set colour
            if ~isempty(options.colour)
                try
                colourvals = GetTubeValues(obj, options.colour, options.colouridx);
                [cdata, rgb] = ColourIndex(obj, colourvals, name=options.colour);                
                catch
                    warning(['Attempted to colour image. Failed ' ...
                        'due to incomplete region definition for all tubes.'])
                    rgb = [];

                end
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

            % plot each color vol as isosurface individually

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
                'FaceColor', 'k',...
                'EdgeColor', 'none');
                hold on
            end

            obj.vol3daxes()
            hold off
        end

        function PlotSpline(obj, options)
            % Plot the splines of all :attr:`tubes` objects.
            %
            % .. todo: 
            %   * stub
            %   * consider weight parameter.
            %
            % Args:
            %   gen(int): *OPTIONAL* `default =
            %       max([obj.tubes.generation])` plot up to which generation.
            %       All by default.
            %   colour(char) = *OPTIONAL* `default = 1`. Set variable for
            %       colour labelling.
            %       if `char` must be an :class:`tube` property, :class:`tube` `stats` or `region` field.
            %       e.g. `'generation'`, `'arclength'`, `'lobe'`.
            %       Can also be vector of length equal to number of tubes
            %       in order. 
            %   colouridx(scalar) = *OPTIONAL* `default = 1` Index of
            %       chosen property in `colour`.    
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.PlotSpline();
            %
            % .. |network_Network_PlotSpline| image:: figs/network_plotspline.png
            %    :width: 400
            %    :alt: figure plot - Network PlotSpline
            %
            % |network_Network_PlotSpline|
            %
            arguments
                obj
                options.gen = max([obj.tubes.generation])
                options.colour = ''
                options.colouridx = 1
            end

            % set up reduced graph
            vis_Glink_logical = [obj.tubes.generation] <= options.gen;
            vis_tubes = obj.tubes(vis_Glink_logical);

            % set colour
            if ~isempty(options.colour)
                try
                colourvals = GetTubeValues(obj, options.colour, options.colouridx);
                [cdata, rgb] = ColourIndex(obj, colourvals, name=options.colour);
                hold on
                catch
                    warning(['Attempted to colour image. Failed ' ...
                        'due to incomplete region definition for all tubes.'])
                end
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

        function s = OrthoView(obj, type, options)
            % View volume using MATLAB's inbuilt othogonal viewer.
            %
            % Call MATLAB's orthosliceViewer for a property that is a 3D
            % array of this class. All params are inferred from the object
            % properties.
            %
            % .. todo::
            %   * look at exposing display range to user.
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
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.OrthoView();
            %
            % .. |network_Network_OrthoView| image:: figs/network_orthoview.png
            %    :width: 400
            %    :alt: figure plot - Network OrthoView
            %
            % |network_Network_OrthoView|
            %

            arguments
                obj
                type {mustBeMember(type,{'source','seg'})} = 'source'
                options.displayrange = [-inf inf]
            end

            % get volume
            volout = ParseVolOut(obj, type=type);
            volout = permute(volout, [2,1,3]);
            volout = flip(volout, 3);

            if islogical(volout)
                volout = single(volout);
            end

            automin = min(volout(:));
            automax = max(volout(:));

            % display range
            switch type
                case 'source'
                    if options.displayrange(1) == -inf
                        options.displayrange(1) = automin;
                    end
                    if options.displayrange(2) == inf
                        options.displayrange(2) = automax;
                    end
                case 'seg'
                    options.displayrange = [automin, automax];
            end


            % display with orthoview

            s = orthosliceViewer(volout, 'DisplayRange', options.displayrange,...
                'DisplayRangeInteraction','off', ...
                'ScaleFactors',obj.voxdim, 'CrosshairLineWidth', 0.3);

        end
        
        function [h, out] = Histogram(obj, options)
            % Plot histogram of :attr:`tubes`.
            %
            %
            % .. note::
            %
            %   Consider case where label contains NaNs.
            %
            % .. warning::
            %
            %   Behaviour untested for label containing NaNs.
            %
            %
            % Args:
            %   label = *OPTIONAL* `default = 'generation'` set variable to count.
            %       if `char` must be an :class:`tube` property, :class:`tube` `stats` or `region` field.
            %       e.g. `'generation'`, `'arclength'`, `'lobe'`.
            %       Can also be vector of length equal to number of tubes
            %       in order. 
            %   labelidx(scalar) = *OPTIONAL* `default = 1`. Index of
            %       chosen property in `label`.
            %   print(bool) = *OPTIONAL* `default = false`. print frequency
            %       table.
            %  
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.PlotGen();
            %
            %
            arguments
                obj
                options.label = 'generation'
                options.labelidx = 1
                options.print = false
            end
            
            % get index for each branch
            outlabel = GetTubeValues(obj, options.label, options.labelidx);
            outlabel = outlabel';

            % get unique values
            outlabel_unique = unique(outlabel);

            % set up cell for each unique val
            rownames = compose('%d', outlabel_unique);
            rownames = cellstr(string(rownames));

            % count for all gens
            [count, ~] = histcounts(outlabel, [outlabel_unique; outlabel_unique(end)+1]);
            % show figure result as bar chart
            h = bar(outlabel_unique, count');
            title(['Number of tubes per ', options.label])
            xlabel(options.label)
            ylabel('count')
            grid on

            if nargout > 1 || options.print == true
                reporttable = table(count', 'RowNames', rownames, 'VariableNames', {options.label});
            end

            if options.print == true
                disp(['Number of tubes per ', options.label])
                disp(reporttable)
            end

            if nargout > 1
                out = reporttable;
            end
            
        end
        % Data IO
        
        function ExportOrthoPatches(obj, path, casename, options)
            % export perpendicular slice patches of all :attr:`tubes`.
            %
            % See :meth:`tube.Tube.ExportOrthoPatches` for more detailed
            % information.
            %
            % Args:
            %   path(str): path to directory to save the exported patches.
            %       The directory will be created if it doesn't already 
            %       exist.
            %   casename(char): casename to use as prefix to each patch
            %       slice name.
            %   genrange(int): *OPTIONAL* `default = [0 inf]` interval
            %   range of generations that should be considered.
            % 
            % Example:
            %   >>> run CA_base.m;
            %   >>> AQnet.ExportOrthoPatches('patches','example')
            %

            arguments
                obj
                path
                casename
                options.genrange = [0 inf]
            end
                
            assert(all(size(options.genrange) == [1 2]), ['genrange ' ...
                'must be of size [1,2]'])
            % set upper level to max possible if set to inf
            if options.genrange(2) == inf
                options.genrange(2) = max([obj.tubes.generation]);
            end
            
            % get list of tubes that meet criteria
            tubelistidx = [obj.tubes.generation] >= options.genrange(1) & ...
                [obj.tubes.generation] <= options.genrange(2);
            tubelist = obj.tubes(tubelistidx);

            for ii = progress(1:length(tubelist), 'Title', 'ExportOrthoPatches')
                tubelist(ii).ExportOrthoPatches(path, casename);
            end
        end

        function obj = ExportCSV(obj, path)
            % Export all tubes and their holistic properties and measurements to csv file.
            %
            % A list of properties from each tube :attr:`ID`, :attr:`parent`, 
            % :attr:`children` and :attr:`generation` and :attr:`method`. 
            % Stats based on measurements from :attr:stats and :attr:region saved into a single 
            % csv file where each tube is represented by a row. 
            %
            % .. todo::
            %   * showcase example
            %
            % Args:
            %   path(char): filename to save as.
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.ExportCSV(example_tubenetwork);
            %
            %

            % parse path to end in csv
            path = parse_filename_extension(path, '.csv');

            % properties to add
            tubeprops = ["ID", "parent", "children", "generation", "method"];
            
            % stuct properties to add
            structprops = ["stats", "region"];
            
            % make new struct per tube and add to table
            exporttable = table;

            for ii = 1:length(obj.tubes)
                atube = obj.tubes(ii);
                % add named tube properties
                rowstruct = struct;
                for prop = tubeprops
                    theprop = atube.(prop);
                    if isa(theprop,'Tube')
                        rowstruct.(prop) = [theprop.ID];
                    else
                        rowstruct.(prop) = theprop;
                    end

                end

                % add tube subproperties
                for prop = structprops
                    subprops = fieldnames(atube.(prop));
                    for sprop = string(subprops)'
                        rowstruct.(strcat(prop,'_',sprop)) = atube.(prop).(sprop);
                    end
                end
                
                % add new row
                row = struct2table(rowstruct, 'AsArray',true);
                exporttable = AQ_vertcat(exporttable, row);
            end

            % write to csv
            writetable(exporttable, path)

        end
    end

    methods(Static)
        function obj = load(path)
            % internal method called when tubenetwork .mat file is loaded.
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
