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
    %   header(struct): header information from source
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
        header
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
        isloops
    end

    methods
        % init
        function obj = TubeNetwork(skel, options)
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
            %   skel (3darray): Binary airway centreline in the
            %       same grid space as CT. Dimensions must match with CT
            %   source (3darray): *OPTIONAL* CT loaded from nifti using niftiread.
            %   header (struct): *OPTIONAL* CT metadata loaded from nifti using
            %       niftiinfo.
            %   seg (3darray): *OPTIONAL* Binary airway segmentation in the
            %       same grid space as CT. Dimensions must match with CT.
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
            skel (:,:,:) logical
            options.source (:,:,:)
            options.header struct
            options.seg (:,:,:) logical
            options.fillholes logical = 1
            options.largestCC logical = 0
            options.spline_sample_sz = nan
            options.plane_sample_sz = nan
            options.max_plane_sz = nan
            options.reorient logical = 1
            options.voxdim (1,3)
            options.originmethod = 'largest';
            end

            assert(all(~skel,'all') == false, 'Skel appears to be empty.')

            if isfield(options,'source')
            assert(all(size(options.source)==size(skel)),['Size of source ',num2str(size(options.source)), ...
                ' differs from skel ', num2str(size(skel))])
            end

            if isfield(options,'seg')
            assert(all(size(options.seg)==size(skel)),['Size of seg ', ...
                size(options.seg),' differs from skel ', num2str(size(skel))])
            end

            if isfield(options,'header')
                obj.header = options.header;
            end
            
            % process skeleton
            robustskel = ParseSeg(skel, fillholes=false, ...
                largestCC=options.largestCC);

            % process segmentation
            if isfield(options,'seg')
                robustseg = ParseSeg(options.seg, fillholes=options.fillholes, ...
                    largestCC=options.largestCC);
            end

            % reorient volumes and get properties
            if options.reorient == true && isfield(options,'header')
                [obj.skel, obj.voxdim] = ReorientVolume(robustskel, obj.header);
                if isfield(options,'source')
                    obj.source = ReorientVolume(options.source, obj.header);
                end
                if isfield(options,'seg')
                    obj.seg = ReorientVolume(robustseg, obj.header);
                end
            else
                obj.skel = robustskel;
                if isfield(options,'seg')
                    obj.seg = robustseg;
                end
                if isfield(options,'source')
                    obj.source = options.source;
                end
            end
            
            % use provided voxdim if exists
            if isfield(options,'voxdim')
                obj.voxdim = options.voxdim;
            elseif ~isfield(options,'header')
                obj.voxdim = [1, 1, 1];
            end

            disp('Dimensions in physical units, usually millimetres.')
            disp(['Voxel dimensions = ', ...
                num2str(obj.voxdim)])

            % identify cropped size by seg if available
            if isfield(options,'seg')
                obj.lims = CropVol(obj.seg);
            else
                obj.lims = CropVol(obj.skel);
            end

            % crop all images
            if isfield(options,'seg')
                obj.seg = (CropVol(obj.seg, obj.lims));
            end
            if isfield(options,'source')
                obj.source = CropVol(obj.source, obj.lims);
            end
            obj.skel = CropVol(obj.skel, obj.lims);

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
            [~, glink, ~] = Skel2Digraph(obj, 'topnode');

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
            
            % classify segmentation to tubes
            if isfield(options,'seg')
                classedseg = obj.ClassifySegmentationTubes();
            else
                classedseg = obj.ClassifySkeletonisationTubes();
            end

            % set origin node by method
            obj.SetOrigin(options.originmethod,classedseg)

            obj.RunAllTubes('SetGeneration');

            % Compute angles of tubes
            obj.RunAllTubes('ComputeDirections');
            obj.RunAllTubes('ComputeChangeAngle');
            obj.RunAllTubes('ComputeParentAngle');
            obj.RunAllTubes('ComputeSiblingAngle',0);

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
                options.gpu logical = false
            end

            for ii = progress(1:length(obj.tubes), 'Title', 'Making tube patches')

                if strcmp(options.type,'source') || strcmp(options.type,'both')
                    assert(~isempty(obj.source),['No source volume to derive patches.' ...
                        ' May need to change `type`, currently set to: ', options.type, '.'])
                    MakePatchSlices(obj.tubes(ii), obj.source, ...
                        type='source', usesegcrop=options.usesegcrop,...
                        method=options.method, gpu=options.gpu);
                end

                if strcmp(options.type,'seg') || strcmp(options.type,'both')
                    assert(~isempty(obj.seg),['No source volume to derive patches.' ...
                        ' May need to change `type`, currently set to: ', options.type, '.'])                    
                    MakePatchSlices(obj.tubes(ii), obj.seg, type='seg', ...
                        usesegcrop=options.usesegcrop, method=options.method, ...
                        gpu=options.gpu);
                end

            end

        end
        
        function classedskel = ClassifySkeletonisationTubes(obj)
            % create classified skeleton
            classedskel = zeros(size(obj.skel));
            for ii = 1:length(obj.tubes)
                classedskel(obj.tubes(ii).skelpoints) = obj.tubes(ii).ID;
            end
        end

        function classedseg = ClassifySegmentationTubes(obj)
            % Classify :attr:`seg` into its tube components.
            %
            % Classifies every foreground point of :attr:`seg` to the
            % nearest point in :attr:`skel`. Assigning it to the respective
            % :class:`tube.Tube` object in :attr:`tubes`.
            %
        
            % create classified skeleton
            classedskel = obj.ClassifySkeletonisationTubes();

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
            [g, glink, gnode, obj.isloops] = skel_2_digraph(obj.skel, method);
        end

        function ge = TubesAsEdges(obj)
            % makes the digraph based on tube relationships

            % init edgetable
            gn = TubesAsNodes(obj);

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
            
            %%% identify origin node of all graphs
            origintubes = [];
            % loop through every tube
            for ii = 1:length(obj.tubes)
                currenttube = obj.tubes(ii);
                % for current tube, find the origin tube (tube with no parent)
                while ~isempty(currenttube.parent)
                    currenttube = currenttube.parent;
                end
                if ~any(origintubes == currenttube)
                    origintubes = [origintubes, currenttube];
                end
            end

            g = digraph();
            
            % DFS from each zero gen tube, creating edgetable.
            for origintube = origintubes
                % case of isolated tube
                if isempty(origintube.parent) && isempty(origintube.children)
                    g = addnode(g,num2str(origintube.ID));
                    continue
                end
                % independent DFS search from each origin tube.
                edgestosearch = origintube;
                while ~isempty(edgestosearch)
                    current_tube = edgestosearch(1);
                    edgestosearch(1) = [];
                    if ~isempty(current_tube.children)
                        for childtube = current_tube.children
                            g = addedge(g, num2str(current_tube.ID), num2str(childtube.ID));
                            edgestosearch = [edgestosearch, childtube];
                        end
                    end
                end
            end
            % remove node names
            g.Nodes.Name = '';
        end
    
        function obj = SetOrigin(obj, method, classedseg)
        % set the origin node by desired method
        if isnumeric(method)
            % numerical origin in RAS provided.
            if length(method) == 3 && numel(method) == 3
                if isrow(method)
                    % ensure it is a vertical vector
                    method = method';
                end
                % convert RAS origin to IJK in context of image
                ijk = round(obj.RAS2IJK(method));
                if any(ijk < 0) || any(ijk' > size(obj.source))
                    % ensure ijk is within ROI
                    error(['Provided origin is outside of ROI. Got ijk: ', ijk])
                end
                % identify which tube this coordinate is within by
                % seg/skel.
                I = classedseg(ijk(1),ijk(2),ijk(3));
                if I == 0
                    error('origin must be exactly within segmentation (if available, otherwise skeleton)')
                end
            else
                error(['numerical origin must be a vector of length 3. Got ', method])
            end
        elseif strcmp(method, 'topnode') || strcmp(method, 'none')
            return
        elseif strcmp(method, 'largest') && ~isempty(obj.seg)
            % choose tube with most segmentation points if available
            [~,I] = max(cellfun(@length,{obj.tubes.segpoints}));
        elseif strcmp(method, 'largest') && isempty(obj.seg)
            % choose tube with longest skeleton
            [~,I] = max(cellfun(@length,{obj.tubes.skelpoints}));
        else
            error(['Set valid origin method. Got ', method])
        end
        obj.tubes(I).SetRoot()

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
                try
                obj.tubes(ii).(tubefunc)(varargin{:});
                catch e
                    warning([char(tubefunc), ' failed for tube index ', num2str(ii), ' .'])
                    fprintf(2,'error identifier :\n%s\n',e.identifier);
                    fprintf(2,'There was an error! The message was:\n%s\n',e.message);
                end
            end
        end

        function I = S2I(obj,I1,I2,I3)
            % Fast specific implementation of MATLAB's `sub2ind
            % <https://uk.mathworks.com/help/matlab/ref/sub2ind.html>`_.
            %

            I = S2I3(size(obj.skel),I1,I2,I3);
        end

        function [I1,I2,I3] = I2S(obj,I)
            % Fast specific implementation of MATLAB's `ind2sub
            % <https://uk.mathworks.com/help/matlab/ref/ind2sub.html>`_.
            %

            [I1, I2, I3] = I2S3(size(obj.skel),I);
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
                try
                    label_all = [eval(object_eval).(labelname)];
                catch % incase tube doesn't have the label
                    label_all = NaN;
                end
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

        function [statsopt, regionopt] = VisualisationOptions(obj)
            % get visualisation options from :attr:`tubes` by :attr:`tube.Tube.stats` and :attr:`tube.Tube.region`.
            statsopt = [];
            regionopt = [];

            % get stats options
            for ii = 1:length(obj.tubes)
                if ~isempty(obj.tubes(ii).stats)
                    current_stats = fieldnames(obj.tubes(ii).stats);
                    statsopt = [statsopt, setdiff(current_stats,statsopt)];
                end
                if ~isempty(obj.tubes(ii).region)
                    current_region = fieldnames(obj.tubes(ii).region);
                    regionopt = [regionopt, setdiff(current_region,regionopt)];
                end
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
                options.colorbar = true
            end

            assert(length(obj.tubes) == length(vals), ['Should have same number ' ...
                'of vals input as there are tubes. Got ', num2str(length(vals)), ...
                ', but expected ', num2str(length(obj.tubes))])
            
            % convert labels into numbers
            if isempty(options.types)
                if isfield(obj.regioncategories,options.name)
                    types = obj.regioncategories.(options.name);
                else
                    types = unique(vals);
                end
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
            if options.colorbar
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

        function [h, ge] = Plot(obj, options)
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
            %   weightfactor(float) = *OPTIONAL* `default = NaN`
            %     determines the highest scaling of the linethickness. if
            %     Nan then no scaling is applied.
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
                options.weightfactor = NaN
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
            if ~isnan(options.weightfactor)
                scale = options.weightfactor/max_thick;
            else
                scale = 1;
            end

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
            %   `'seg'` or `'skel'`. which to plot surface.
            %  smooth_sz(int): *OPTIONAL* `default = 0` size of gaussian 
            %   smoothing kernel. 0 for no smoothing.
            %
            % .. note:
            %   Plotting multiple patches (>20) significantly reduces
            %   performance. This function attempts to group
            %   segments of the same colour to plot together.
            % 
            % .. warning:
            %   Choosing a colour option with several colours will cause this
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
                options.smooth_sz = 0
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
                % smoothing
                if options.smooth_sz > 0
                    U = smooth3(U,'gaussian',options.smooth_sz);
                end
                patch(isosurface(U),...
                'FaceAlpha', options.alpha,...
                'FaceColor', rgb(ii,:),...
                'EdgeColor', 'none');
                hold on
                end
            else
                if options.smooth_sz > 0
                    V = smooth3(V,'gaussian',options.smooth_sz);
                end
                patch(isosurface(V),...
                'FaceAlpha', options.alpha,...
                'FaceColor', 'c',...
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
        
        function [h, out] = Histogram2(obj, options)
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
            % .. todo::
            %
            %   Remove redundant variables.
            %
            % Args:
            %   label = *OPTIONAL* `default = 'generation'` set variable to count.
            %       if `char` must be an :class:`tube` property, :class:`tube` `stats` or `region` field.
            %       e.g. `'generation'`, `'arclength'`, `'lobe'`.
            %       Can also be vector of length equal to number of tubes
            %       in order. 
            %   labelidx(scalar) = *OPTIONAL* `default = 1`. Index of
            %       chosen property in `label`.
            %   print(bool) = *OPTIONAL* `default = false`. Print frequency
            %       table.
            %   exact(bool) = *OPTIONAL* `default = false`. Bin values by
            %       their exact number if label is numerical.
            %  
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.Histogram();
            %
            %
            arguments
                obj
                options.label = 'generation'
                options.labelidx = 1
                options.print = false
                options.exact = false
            end
            
            % get index for each branch
            outlabel = GetTubeValues(obj, options.label, options.labelidx);
            outlabel = outlabel';

            % get unique values
            outlabel_unique = unique(outlabel);
            categoricalinput = 0;
            % set up cell for each unique val
            if isnumeric(outlabel_unique)
                rownames = compose('%d', outlabel_unique);
                rownames = cellstr(string(rownames));
                % 2nd variable defines the right side of each bin only
            else
                % must be categorical, convert
                categoricalinput = 1;
                rownames = outlabel_unique;
                outlabel = categorical(outlabel);
                outlabel_unique = categorical(outlabel_unique);
            end

            if ~options.exact || categoricalinput
                % determine counts by binning algorithm if numeric
                [count, edges] = histcounts(outlabel);
            else
                % use exact bins
                [count, edges] = histcounts(outlabel, [outlabel_unique; outlabel_unique(end)+1]);
            end
            
            if categoricalinput
                % convert to categoric
                edges = categorical(edges);
            else
                % remove right most edge if numerical.
                edges = edges(1:end-1);
            end

            % count for all vals
            % show figure result as bar chart
            h = bar(edges, count');
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
 
        function h = Histogram(obj, options)
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
            % .. todo::
            %
            %   Remove redundant variables.
            %
            % Args:
            %   label = *OPTIONAL* `default = 'generation'` set variable to count.
            %       if `char` must be an :class:`tube` property, :class:`tube` `stats` or `region` field.
            %       e.g. `'generation'`, `'arclength'`, `'lobe'`.
            %       Can also be vector of length equal to number of tubes
            %       in order. 
            %   labelidx(scalar) = *OPTIONAL* `default = 1`. Index of
            %       chosen property in `label`.
            %   region(char) = *OPTIONAL* `default = ''`. Region to
            %       group label values by.
            %  
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.Histogram2();
            %
            %
            arguments
                obj
                options.label = 'generation'
                options.labelidx = 1
                options.region = ''
            end
            
            % get index for each branch
            outlabel = GetTubeValues(obj, options.label, options.labelidx);
            outlabel = outlabel';
            
            if ~isempty(options.region)
                regions = GetTubeValues(obj, options.region, 1);
            else
                regions = cell(1, length(outlabel));
                regions(:) = {'all'};
            end
            
            if iscell(regions)
                % make any possible nans into chars
                [regions{~cellfun(@ischar,regions)}] = deal('none');
                % convert to categorical for histogram function
                regions = categorical(regions);
            end
            regions_unique = unique(regions);

            % count for all vals
            % show figure result as bar chart
            tcl = tiledlayout('flow');
            for iregion = regions_unique
                nexttile
                h = histogram(outlabel(regions == iregion));
                title(iregion)
                xlabel(options.label)
                ylabel('count')
            end
            title(tcl,options.label)
            grid on


        end
        % Data IO

        function [ijk] = RAS2IJK(obj,ras)
            % RAS2IJK converts from RAS to IJK coordinates and crops.
            % see :func:`RAS2IJK` for more details.
            %
            
            % apply affine transform
            ijk_whole = ras_2_ijk(ras, obj.header);

            % crop
            ijk = ijk_whole - obj.lims(:,1) + 1;

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
            
        % export to csv
        writetable(node_table, path)

        end
        
        function ExportSlicerLine(obj, path, options)
            arguments
                obj
                path
                options.n = length(obj.tubes);
                options.colour = ''; % turquoise
                options.colouridx = 1
            end

            % Export tubes as JSON line markups for viewing in slicer
            
            % ensure correct extension on path
            path = parse_filename_extension(path, '.mrk.json');
            
            % attempt to set colour
            if ~isempty(options.colour)
                try
                colourvals = GetTubeValues(obj, options.colour, options.colouridx);
                [cdata, rgb] = ColourIndex(obj, colourvals, name=options.colour, colorbar=false);                
                catch
                    warning(['Attempted to colour image. Failed ' ...
                        'due to incomplete region definition for all tubes.'])
                    rgb = [];

                end
            else
                rgb = [];
            end


            % process tubes
            tubemarkups = cell(options.n,1);
            for ii = 1:options.n
                if ~isempty(rgb)
                    tubemarkups{ii,1} = obj.tubes(ii).ExportSlicerLine(rgb(cdata(obj.tubes(ii).ID),:));
                else
                    tubemarkups{ii,1} = obj.tubes(ii).ExportSlicerLine();
                end
            end
            
            % make json sting
            jsondict = struct();
            jsondict.schema = "https://raw.githubusercontent.com/slicer/slicer/master/Modules/Loadable/Markups/Resources/Schema/markups-schema-v1.0.3.json#";
            jsondict.markups = tubemarkups;
            jsonstr = jsonencode(jsondict,PrettyPrint=true);
            jsonstr = replace(jsonstr,'"schema":','"@schema":');
            
            % save
            fid = fopen(path,'w');
            fprintf(fid,'%s',jsonstr);
            fclose(fid);
        end

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
        
        function grid_preview_measures(obj, n_rows, n_cols)
            arguments
                obj
                n_rows
                n_cols
            end
            % export grid of preview of patches, with measures
            %

            % map patch count to continuous index
            for id = 1:length(obj.tubes)
                tube = obj.tubes(id);
                n_tube_patches = size(tube.source,3);
                map = [ones(n_tube_patches,1)*id,(1:n_tube_patches)'];
                if id == 1
                    mapping = map;
                else
                    mapping = [mapping; map];
                end
            end
            npatches = length(mapping);
            % sample n desired patches
            randi = randperm(npatches);
            
            f = gcf;
            % force figure to be square
            f.Position(3:4) = [800, 800];
            tiledlayout(n_rows,n_cols,"TileSpacing","none")
            
            for k = randi(1:n_rows*n_cols)
                patch_map = mapping(randi(k),:);
                tube = obj.tubes(patch_map(1));
                nexttile;
                tube.View(patch_map(2));
                axis off
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
