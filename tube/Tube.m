% Lead Author: Ashkan Pakzad 2022. ashkanpakzad.github.io.
% See https://github.com/ashkanpakzad/AirQuant for more information.

classdef Tube < matlab.mixin.SetGet
    % Initialise the Tube class object.
    %
    % description
    %
    % .. todo: add documentation to this function
    %
    % Args:
    %   network (:class:`AirQuant.network`): tube network object that this
    %     tube is a subset
    %   relatives (struct): related tube objects e.g. `parent` to  tube
    %     objects
    %   skelpoints (int, vector): list of linear indexed points that make
    %     up the tube's centreline
    %   spline (struct): polynomial that describes the spline
    %   source (float, array): tube patches of source image along
    %   spline at interval specifed by
    %     :attr:`patchprop.arclength`
    %   seg (float, array): tube patches of seg image along spline
    %     at interval specifed by :attr:`patchprop.arclength`
    %   patchprop (struct): property per tube patch slice given as list
    %     e.g :attr:`patchprop.arclength`
    %   stats (struct): stats property for current tube
    %
    %
    % see also

    properties
        parent = []
        children = []
        generation
        network
        relatives
        skelpoints
        spline
        source
        seg
        patchprop
        prunelength = []
        stats
        region
        savename
    end
    properties (SetAccess = protected)
        ID
        measures = []
        method
        diameters = []
        areas = []
        volumes = []
    end

    methods
        function obj = Tube(network, skelpoints, ID)
            % tube representing anatomical structure
            %
            % .. todo:
            %   * make tubes saveable and loadable and that these operations
            %   can be done independently of the network object.
            %   * to make independent, it may be necessary to remove the
            %   network as a property. it will be nessesary to save the
            %   keep the source size.
            %
            % Args:
            %   network : network object to which this tube belongs.
            %   skelpoints : list of points that make up the tube's
            %   skeleton to the network's source image.
            %   localsegpoints : list of points that make up the tube's
            %   segmentation to the network's source image.
            %   parent(cell): `OPTIONAL` parent tube(s) that directs to this tube
            %       object.
            %   child(cell): `OPTIONAL` child tube(s) that directs from this
            %       tube object.
            %

            arguments
                network
                skelpoints (1,:)
                ID (1,1)
            end

            % definitions
            obj.network = network;
            obj.skelpoints = skelpoints;
            obj.ID = ID;

            obj.relatives = struct;
            obj.patchprop = struct;
            obj.stats = struct;
            obj.region = struct;

            obj.MakeSpline();
            obj.ComputeSplinePoints();

        end

        function obj = SetChildren(obj, tube)
            % Set relative to current tube object.
            %
            % desc
            % .. todo: add documentation to this function
            %
            % Args:
            %   relativetube (:class:`tube`): the tube to set
            %       relation to.
            %   relation (string): relation name. common
            %       "parent" or "child".
            %

            obj.children = [obj.children tube];
            tube.SetParent(obj)
        end

        function obj = SetParent(obj, tube)
            % Set relative to current tube object.
            %
            % desc
            % .. todo:
            %   * add documentation to this function
            %   * set children tube of parent without being stuck in loop.
            %
            % Args:
            %   relativetube (:class:`tube`): the tube to set
            %       relation to.
            %   relation (string): relation name. common
            %       "parent" or "child".
            %
            obj.parent = [obj.parent tube];
        end

        function obj = SetGeneration(obj)
            % Set relative to current tube object.
            %
            %
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   relativetube (:class:`tube`): the tube to set
            %       relation to.
            %   relation (string): relation name. common
            %       "parent" or "child".
            %
            %
            currentbranch = obj;
            count = 0;
            while ~isempty(currentbranch.parent)
                count = count + 1;
                currentbranch = currentbranch.parent;
                if length(currentbranch.parent) > 1
                    obj.generation = NaN;
                    warning('Multiple parents, not possible to set generation')
                    return
                end
            end
            obj.generation = count;

        end

        function obj = SetRegion(obj, regiontype, value)
            % Set region classifcation of tube.
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   relativetube (:class:`tube`): the tube to set
            %       relation to.
            %   relation (string): relation name. common
            %       "parent" or "child".
            %
            %
            obj.region = setfield(obj.region, regiontype, value);
        end

        % spline related
        function obj = MakeSpline(obj, options)
            % fits a spline to the centreline of the tube.
            %
            % Using the object property :attr:`skelpoints`, a polynomial
            % spline is fit to this list of points. A moving average is
            % used to smooth the skeletal points, it can also use the
            % parent tube to initialise this moving average. For a better
            % understanding of the spline output see `cscvn`.
            % Based on original function by Kin Quan 2018
            %
            % .. warning:
            %   the order of skeletal points affects outcome e.g. reversing
            %   the order of the skeleton points would reverse the spline
            %   direction.
            %
            % Args:
            %   useparent(bool): `optional` uses the parent tube skel
            %       points if available to initialise the moving average.
            %       True by default.
            %
            % Return:
            %   obj: object :attr:`spline` property updated. A `cscvn`
            %       struct object.
            %
            arguments
                obj
                options.useparent logical = true
            end

            % get linear indexed points of previous branch if available.
            if options.useparent == true && isfield(obj.relatives,'parent')
                parent_points = obj.relatives.parent.skelpoints;
                [x_p1, y_p1, z_p1] = ind2sub(size(obj.network.source), parent_points);
            else
                x_p1 = []; y_p1 = []; z_p1 = [];
            end

            % get current tube points
            [x_p2, y_p2, z_p2] = I2S(obj, obj.skelpoints);
            x_point = [x_p1, x_p2];
            y_point = [y_p1, y_p2];
            z_point = [z_p1, z_p2];

            %Smooth all points using moving average
            voxel_sz = obj.network.voxdim;
            smooth_x = smooth(x_point*voxel_sz(1),11, 'moving');
            smooth_y = smooth(y_point*voxel_sz(2),11, 'moving');
            smooth_z = smooth(z_point*voxel_sz(3),11, 'moving');

            % extract just current smoothed points
            csmooth_x = smooth_x(length(x_p1)+1:end);
            csmooth_y = smooth_y(length(x_p1)+1:end);
            csmooth_z = smooth_z(length(x_p1)+1:end);

            %Complete smooth data
            smooth_data_points = [csmooth_x csmooth_y csmooth_z]';

            %Generating the spline
            obj.spline = cscvn(smooth_data_points);
        end

        function obj = ComputeSplinePoints(obj, options)
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
            arguments
                obj
                options.sample_interval (1,1) = obj.network.spline_sample_sz
            end
            assert(~isempty(obj.spline), 'spline is empty, see method MakeSpline')

            % get spline points by set interval
            [obj.stats.arclength, obj.patchprop.parapoints, ...
                obj.patchprop.arcpoints] = ...
                Compute_Spline_Points(obj.spline, options.sample_interval);

            % save stats measurement using derived spline points.
            obj = ComputeTortuosity(obj);
        end

        % perpendicular slice interpolation
        function obj = MakePatchSlices(obj, vol, options)
            % Constructs perpendicular images as if travelling along an
            % airway segment in CT image and Segmentation.
            % short desc
            %
            % long desc
            %
            % .. todo:
            %   * add documentation to this function
            %   * Heavily reliant on the network class structure.
            %       Consider decoupling network in this function.
            %   * Consider incorporating matlabs obliqueslice function.
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            arguments
                obj
                vol (:,:,:)
                options.type char {mustBeMember(options.type,{'source','seg'})} = 'infer'
                options.usesegcrop logical = false
            end

            assert(~isempty(obj.spline), 'spline is empty, see method MakeSpline')
            if nargin < 3

            end

            warn_arg_comb = 'invalid vol and type arg combination see documentation';

            if isa(vol,'numeric')
                assert(~strcmp(options.type,'seg'), warn_arg_comb)
                if ~strcmp(options.type,'source')
                    warning('Assumed source image, saving to reformed source.')
                end
                reformedproperty = 'source';
            end

            if isa(vol,'logical')
                assert(~strcmp(options.type,'source'), warn_arg_comb)
                if ~strcmp(options.type,'seg')
                    warning('assumed seg image, saving to reformed seg.')
                end
                reformedproperty = 'seg';
            end
            
            % cast volume to single for interpolation
            vol = single(vol);

            if options.usesegcrop == true
                obj.patchprop.approx_diameter = NaN(size(obj.patchprop.arclength));
            end

            % set up slice store
            reformedimages = cell(length(obj.patchprop.parapoints),1);
            for i = 1:length(obj.patchprop.parapoints)
                % Compute Normal Vector per spline point
                [normvec, point] = spline_normal(obj.spline, ...
                    obj.patchprop.parapoints(i));

                % plane size
                max_sz = obj.network.max_plane_sz;
                plane_sz = max_sz+1;
                if options.usesegcrop == true
                    assert(obj.network.plane_scaling_sz > 0, ...
                        'obj.network.plane_scaling_sz must be real positive')
                    scaling_sz = obj.network.plane_scaling_sz;
                    % get approx size from distance map of seg
                    obj.patchprop.seg_diameter(i) = ApproxSegDiameter(obj, point, i);
                    plane_sz = ceil(approx_diameter*scaling_sz);

                end
                % use max plane size if current plane size exceeds it
                if plane_sz > max_sz
                    plane_sz = max_sz;
                end

                % Interpolate Perpendicular Slice per spline point
                reformedimages{i,1} = PlaneInterpVol(vol, ...
                    obj.network.voxdim, point, normvec, ...
                    plane_sz=plane_sz, ...
                    sample_sz=obj.network.plane_sample_sz, ...
                    offgrid_val=0);
            end

            % Save plane images
            if strcmp(reformedproperty, 'source')
                obj.source = reformedimages;
            elseif strcmp(reformedproperty, 'seg')
                obj.seg = reformedimages;
            end
        end

        function segdiameter = ApproxSegDiameter(obj, sourcepoint)
            % short desc
            %
            % long desc
            %
            % .. todo: * add documentation to this function
            %   * consider diameter conversion to mm
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %

            % Convert CT_point mm back to voxel ind
            vox_point = sourcepoint'./obj.network.voxdim;
            % find nearest skel point to voxpoint
            assert(~isempty(obj.skelpoints),'need skelpoints defined')
            P = obj.skelpoints;
            k = dsearchn(P,vox_point);
            % Get radius and convert to diameter
            segdiameter = obj.network.Dmap(P(k,1),P(k,2),P(k,3))*2;
            % incase of edge case, unit radius
            if segdiameter == 0
                segdiameter = 2;
            end
        end
    
        % stats
        function obj = ComputeEucLength(obj)
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
            assert(isfield(obj.patchprop,'parapoints'), 'no parapoints computed, see method ComputeSplinePoints')
            [~, point_1] = spline_normal(obj.spline, ...
                obj.patchprop.parapoints(1));
            [~, point_end] = spline_normal(obj.spline, ...
                obj.patchprop.parapoints(end));
            obj.stats.euclength = norm(point_end - point_1);
        end

        function obj = ComputeArcLength(obj)
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
            obj.stats.arclength = Compute_Spline_Points(obj.spline);
        end

        function obj = ComputeTortuosity(obj)
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
            if ~isfield(obj.stats,'arclength')
                obj = ComputeArcLength(obj);
            end
            obj = ComputeEucLength(obj);
            % arclength / euclidean length
            obj.stats.tortuosity = obj.stats.arclength./obj.stats.euclength;
            assert(obj.stats.tortuosity >= 1, 'Impossible to get a tortuosity > 1')
        end
        
        function meanval = ComputeMean(obj, patchprop, trim)
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
            if nargin < 3
                trim = 0;
            end
            assert(isfield(obj.patchprop, patchprop), 'input must be named variable in obj.patchprop.')
            
            % prune the two variables
            var = obj.PruneMeasure(getfield(obj.patchprop,patchprop));

            % compute average
            meanval = trimmean(var, trim);
            obj.stats.('patchprop').trimmean = meanval;
            obj.stats.('patchprop').trim = trim;
        end

        function intrataperval = ComputeIntrataper(obj, patchprop)
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
            assert(isfield(obj.patchprop, patchprop), 'input must be named variable in obj.patchprop.')

            % prune the two variables
            al = obj.PruneMeasure(obj.patchprop.arclength);
            var = obj.PruneMeasure(getfield(obj.patchprop,patchprop));
            % fit bisquare method
            coeff = robustfit(al, var,'bisquare');
            % compute intra-branch tapering as percentage
            intrataperval = -coeff(2)/coeff(1) * 100;
            obj.stats.(patchprop).intrataper = intrataperval;
        end
        
        function intertaperval = ComputeIntertaper(obj, patchprop, trim)
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
            assert(isfield(obj.patchprop, patchprop), 'input must be named variable in obj.patchprop.')
            assert(length(obj.parents)==1, 'must only have one parent tube.')

            % compute means
            parentmean = obj.parent.ComputeMean(patchprop, trim);
            currentmean = obj.ComputeMean(patchprop, trim);

            % compute interbranch tapering as percentage
            intertaperval = (parentmean - currentmean)/(parentmean) * 100;
            obj.stats.(patchprop).intertaperval = intertaperval;
        end
        
        % measures
        function obj = Measure(obj, classmethod, varargin)
            % Call desired method to make measurement on tube patch slices.
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

            % reset measures property
            obj.measures = [];
            obj.method = classmethod;
            % reset diameter and area
            obj.diameters = [];
            obj.areas = [];
            % call measure method
            varargin = {obj.network.plane_sample_sz, obj, varargin{:}};
            classmeasures = feval(classmethod, varargin{:});
            obj.measures = classmeasures.measures;
            obj.diameters = classmeasures.OutputDiameter();
            obj.areas = classmeasures.OutputArea();
        end

        % visualisation
        function h = PlotMeasure(obj, patchprop, varargin)
            % plot patchprop measure
            %
            % desc
            %
            % .. todo:
            %   * add documentation to this function
            %   * Needs attention
            %
            % Args:
            %  relativetube (:class:`tube`): the tube to set
            %   relation to.
            %  relation (string): relation name. common
            %   "parent" or "child".
            %
            X = obj.patchprop.arclength;
            if isa(patchprop, "char") || isa(patchprop, "string")
                Y = getfield(obj.patchprop, patchprop);
            else
                Y = patchprop;
            end
            assert(size(X) == size(Y), 'patchprop needs to be the same length as number of patch slices.')
            h = plot(X, Y, varargin{:});
        end
        
        function s = View(obj, options)
            % View a series of an airway segment's slices as a volume image
            % stack using MATLAB's inbuilt othogonal 3d viewer.
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

            arguments
                obj
                options.type {mustBeMember(options.type,{'source','seg'})} = 'source'
                options.rings = ones(1,size(obj.measures,1))
                options.ellipses = true
                options.points = false
            end

            assert(size(options.rings,2)==size(obj.measures,1), ...
                ['rings must be same length as number of measurements, got ' ...
                , num2str(size(options.rings,2)), ' instead of ' num2str(size(obj.measures,1),2)])

            % convert from cell stack to 3D array.
            tubearray = tubestack(obj, type=options.type);

            % display with orthoview
            s = orthosliceViewer(tubearray, 'DisplayRangeInteraction','off', ...
                'ScaleFactors',[obj.network.plane_sample_sz, ...
                obj.network.plane_sample_sz, ...
                obj.network.spline_sample_sz], 'CrosshairLineWidth', 0.3);
            % get axial axes
            [ax, ~, ~] = getAxesHandles(s);

            addlistener(s,'CrosshairMoving',@allevents);
            addlistener(s,'CrosshairMoved',@allevents);
            
            % set initial annotation
            obj.UpdateAnnotateOrthoviewer(ax,s.SliceNumbers(3),...
                options.rings,options.ellipses...
                ,options.points)

            function allevents(src,evt)
                evname = evt.EventName;
                switch(evname)
                    case{'CrosshairMoving','CrosshairMoved'}
                        ppos = evt.PreviousPosition(3);
                        cpos = evt.CurrentPosition(3);
                        if ppos ~= cpos
                            obj.UpdateAnnotateOrthoviewer(ax,cpos,...
                                options.rings,options.ellipses...
                                ,options.points)
                        end
                end
            end

        end

        % Data IO
        function SegmentTaperResults = SegmentTaperAll(obj, prunelength)
            % Set relative to current tube object.
            %
            % desc
            %
            % .. todo:
            %   *add documentation to this function
            %   * Needs attention
            %
            % Args:
            %  relativetube (:class:`tube`): the tube to set
            %   relation to.
            %  relation (string): relation name. common
            %   "parent" or "child".
            %
            % high level function to compute the segmental tapering
            % measurement of all airways.
            %

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

        function obj = SaveAwyPatches(obj, prunelength)
            % short desc
            %
            % long desc
            %
            % .. todo:
            %   * add documentation to this function
            %   * Needs attention
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %

            if nargin < 2
                prunelength = [0 0];
            end
            % make directory
            [fPath, saveid, ~] = fileparts(casedir);
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
        
        function toGif(obj, filename, options)
            % View a series of an airway segment's slices as a volume image
            % stack using MATLAB's inbuilt othogonal 3d viewer.
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

            arguments
                obj
                filename
                options.type {mustBeMember(options.type,{'source','seg'})} = 'source'
                options.rings = ones(size(obj.measures,1),1)
                options.ellipses = true
                options.points = false
            end
            
            % parse filename
            filename = parse_filename_extension(filename, '.gif');

            % instantiate orthosliceviewer
            s = obj.View(type=options.type, rings=options.rings, ...
                ellipses=options.ellipses, ...
                points=options.points);
            [ax, ~, ~] = getAxesHandles(s);
            
            set(s,'CrosshairEnable','off');
            sliceNums = 1:length(obj.source);
            for idx = sliceNums
                % Update z slice number and annotation to get XY Slice.
                s.SliceNumbers(3) = idx;
                obj.UpdateAnnotateOrthoviewer(ax,idx,options.rings,...
                    options.ellipses,options.points)

                % Use getframe to capture image.
                I = getframe(ax);
                [indI,cm] = rgb2ind(I.cdata,256);

                % Write frame to the GIF File.
                if idx == 1
                    imwrite(indI,cm,filename,'gif','Loopcount',inf,'DelayTime',0.05);
                else
                    imwrite(indI,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
                end
            end
        end

        function toNii(obj, filename, options)
            arguments
                obj
                filename char
                options.type {mustBeMember(options.type,{'source','seg'})} = 'source'
                options.gz logical = true
            end

            filename = parse_filename_extension(filename, '.nii');

            tubearray = tubestack(obj, type=options.type);
            niftiwrite(tubearray,filename)
            
            if options.gz == true
                gzip(filename)
                delete(filename)
            end
        end
    
        function status = toITKsnap(obj, options)
            % view in itksnap
            %
            % May need to set up enviroments on matlab search path for
            % system terminal. 
            % setenv('PATH', [getenv('PATH') ':/Applications/ITK-SNAP.app/Contents/bin']);
            %
            arguments
                obj
                options.type {mustBeMember(options.type,{'source','seg'})} = 'source'
            end
            
            % set up temp file
            filename = parse_filename_extension(tempname, '.nii');
            obj.toNii(filename, type=options.type, gz=false)
            
            command = ['itksnap ', filename];
            status = system(command);

            if status ~= 0
                error(['Failed to open in ITK-snap, ' ...
                    'please see documentation to check this facility is ' ...
                    'set up correctly.'])
            end

        end
        
        % utilities
        
        function I = S2I(obj,I1,I2,I3)
            % short desc
            %
            % long desc
            %
            % .. todo: add docs
            %
            % Args:
            %   x():
            %
            % Return:
            %   y():
            %

            I = S2I3(size(obj.network.source),I1,I2,I3);
        end

        function [I1,I2,I3] = I2S(obj,I)
            % short desc
            %
            % long desc
            %
            % .. todo: add docs
            %
            % Args:
            %   x():
            %
            % Return:
            %   y():
            %

            [I1, I2, I3] = I2S3(size(obj.network.source),I);
        end
        
        function obj = SetPruneLength(obj, prunelength)
            arguments
                obj
                prunelength (2,1) mustBeNumeric
            end
            obj.prunelength = prunelength;
        end
        
        function prunedprop = PruneMeasure(obj, var)
            if isempty(obj.prunelength)
                obj.prunelength = [0 0];
                warning('Object prunelength not set, using zero pruning.')
            end
            pl = obj.prunelength;
            al = obj.patchprop.arclength;
            assert(length(al) == length(var), 'Input variable must be the same length as arclength.')

            prunebool = (al >= pl & al <= al(end) - pl);
            prunedprop = var(prunebool);
        end

        function tubearray = tubestack(obj,options)
            % short desc
            %
            % long desc
            %
            % .. todo: a
            %   * add documentation to this function
            %   * add version that makes tubestack back to parent
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %

            arguments
                obj
                options.type {mustBeMember(options.type,{'source','seg'})} = 'source'
            end
            % generate an airway's interpolated slices into an array
            % stack.
            tubecell = get(obj,options.type);
            canvas_sz = floor(obj.network.max_plane_sz/obj.network.plane_sample_sz);
            tubearray = zeros([canvas_sz, canvas_sz, length(tubecell)]);
            for slice = 1:length(tubecell)
                image = tubecell{slice,1};
                image_sz = size(image,1);
                min_centre = canvas_sz/2 - image_sz/2;
                max_centre = canvas_sz/2 + image_sz/2;
                tubearray(min_centre+1:max_centre, min_centre+1:max_centre, slice) = image;
            end
        end

    end

    methods (Access = protected)
        function UpdateAnnotateOrthoviewer(obj,ax,pos,rings,showellipses,showpoints)
            % delete old linetype graphics
            axesHandlesToChildObjects = findobj(ax, 'Type', 'Line');
            if ~isempty(axesHandlesToChildObjects)
                delete(axesHandlesToChildObjects);
            end
            for ii = 1:size(obj.measures,1)
                if rings(ii) == 1

                    % get centre displacement
                    canvas_sz = floor(obj.network.max_plane_sz/obj.network.plane_sample_sz);
                    image_sz = size(obj.source{ii,1},1);
                    min_centre = canvas_sz/2 - image_sz/2;
                    obj.measures{ii,pos}.plot(min_centre, ax, showellipses, showpoints);
                end
            end
        end
    end
end
