% Lead Author: Ashkan Pakzad 2022. ashkanpakzad.github.io.
% See https://github.com/ashkanpakzad/AirQuant for more information.

classdef Tube < AirQuant & matlab.mixin.SetGet
    % A tubular structure in a connected tree-like structure.
    %
    % Represents a tube based on some source image for analysis.
    % Where quantitative metrics can be extracted.
    %
    %
    %
    % Args:
    %   source (float): *INHERITED* array of interpolated perpendicular slice patches of source
    %     along the :attr:`tube.Tube.spline`.
    %   seg (float): *INHERITED* array of interpolated perpendicular slice patches of seg
    %     along the :attr:`tube.Tube.spline`.
    %   voxdim (float): *INHERITED* triplet that specifies the dimensions of the
    %     tube source where the third dimension represents the :attr:`tube.Tube.spline` interval.
    %   parent (:class:`tube.Tube`): vector of parent tubes to this tube.
    %   children (:class:`tube.Tube`): vector of child tubes to this tube.
    %   generation (int): scalar, number of ancestor tubes.
    %   network (:class:`network.TubeNetwork`): tube network object that this
    %     tube is a subset.
    %   skelpoints (int): vector of linear indexed points that make
    %     up the tube's centreline
    %   spline (struct): polynomial spline that is created using `cscvn`_.
    %   patchprop (struct): property per tube patch slice given as list
    %     e.g `patchprop.arcpoints`
    %   prunelength (float): two element vector that represents in units of
    %     `patchprop.arclength` how much of the tube to prune at either end.
    %     This is set by :meth:`tube.Tube.SetPruneLength`
    %   stats (struct): stats property for current tube e.g `stats.arclength`
    %   region (struct): region that this tube belongs to e.g. `region.lobe`
    %   ID (int): *PROTECTED* scalar index assigned to tube.
    %   method (char): *PROTECTED* method used to make in-slice measurements e.g. :attr:`tube.Tube.diameters`
    %   diameters (float): *PROTECTED* `m x n` array of diameter measurements made in plane.
    %     where `m` is the number of interpolated slices and `n` is the number of in-plane measurements.
    %   areas (float): *PROTECTED* `m x n` array of area measurements made in plane.
    %     where `m` is the number of interpolated slices and `n` is the number of in-plane measurements.
    %   volumes (float): *PROTECTED* vector of volume measurements, where length is the number of in-plane measurements.
    %
    % .. warning::
    %     Though it is possible to set multiple parents, some functionality will be lost.
    %     e.g. it is not possible to use :meth:`tube.Tube.SetGeneration`
    %
    % .. todo::
    %   * make tubes saveable and loadable and that these operations
    %     can be done independently of the network object.
    %   * decouple from network. i.e. remove network as property, and convert skelpoints to sub index.
    %
    % .. _cscvn: https://uk.mathworks.com/help/curvefit/cscvn.html
    %

    properties
        parent = []
        children = []
        generation
        network
        relatives
        skelpoints
        segpoints
        spline
        patchprop
        prunelength = []
        stats
        region
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
            % Init Tube see :class:`tube.Tube`.
            %
            % Makes the Tube object by defining the skelpoints along the centreline of the tube.
            % This init method automatically interpolates the tube and fits a spline.
            % The network object is required as it specifies certain configuration properties.
            % This function is not intended to be invoked directrly. Only by
            % :meth:`network.TubeNetwork.MakeTubes`.
            %
            % Args:
            %   network : see :attr:`tube.Tube.network`
            %   skelpoints : see :attr:`tube.Tube.skelpoints`
            %   ID : see :attr:`tube.Tube.ID`
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
            obj.voxdim = [obj.network.plane_sample_sz, ...
                obj.network.plane_sample_sz, ...
                obj.network.spline_sample_sz];

            obj.relatives = struct;
            obj.patchprop = struct;
            obj.stats = struct;
            obj.region = struct;

            obj.MakeSpline();
            obj.FindSplinePoints();

        end

        function obj = SetChildren(obj, tube)
            % Set children of tube.
            %
            % This method automatically invokes :meth:`tube.Tube.SetParent`
            % in the children tube specified.
            %
            % Args:
            %   tube (:class:`tube.Tube`): list to set as children.
            %

            obj.children = [obj.children tube];
            tube.SetParent(obj);
        end

        function obj = SetParent(obj, tube)
            % Set parent of tube.
            %
            % .. todo::
            %   * set children tube of parent without being stuck in loop.
            %
            % Args:
            %   tube (:class:`tube.Tube`): list to set as parent.
            %

            obj.parent = [obj.parent tube];
        end

        function obj = SetGeneration(obj)
            % Identify and set the generation of this tube.
            %
            % Counts the number tube anscestors i.e. parents of parents up to
            % parent of generation 0.
            %
            % .. warning::
            %   This will fail if there are more than two parent tubes.
            %
            %
            currentbranch = obj;
            count = 0;
            while ~isempty(currentbranch.parent)
                count = count + 1;
                if currentbranch.generation == 0
                    % effective when the true 0th gen is set for a
                    % particular case.
                    break
                end
                currentbranch = currentbranch.parent;
                if length(currentbranch) > 1
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
            % modifies :attr:`tube.Tube.region` such that
            % tube.Tube.region.regiontype = value.
            %
            %
            % Args:
            %   regiontype (`char`): name of the category to set this region.
            %     e.g. `lobe`
            %   value (): value that this tube should be set for the given
            %     regiontype. e.g. `LLL` for Left lower Lobe.
            %
            %
            obj.region = setfield(obj.region, regiontype, value);
        end

        function obj = SetRegionDescendants(obj, regiontype, value)
            % Set region classifcation of tube for current and all
            % children.
            %
            %
            %
            % Args:
            %   regiontype : See :meth:`tube.Tube.SetRegion`
            %   value : See :meth:`tube.Tube.SetRegion`
            %

            % get all desendants
            for tubeii = obj.Descendants()
                tubeii.SetRegion(regiontype, value);
            end
        end

        function obj = SetRegionGeneration(obj, regiontype)
        % set generation by given region to region property
        %
        % Example:
        %   >>> AQnet.tubes(1).SetRegionGeneration('lobe')
        %   >>> AQnet.tubes(1).region.lobe_gen
        %
        % todo:
        %   * needs attention
        %
        % Args:
        %   regiontype : See :meth:`tube.Tube.SetRegion`
        %
        genname = [regiontype, '_gen'];
        currentube = obj;
        count = 0;
        while ~strcmp(currentube.region.(regiontype), ...
                currentube.parent.region.(regiontype))||...
                isempty(currentube.parent)
            if length(currentube.parent) > 1
                obj.region.(genname) = NaN;
                warning('Multiple parents, not possible to set generation')
                return
            end
            count = count + 1;
            currentube = currentube.parent;
        end
        obj.region.(genname) = count;
        end

        function descendants = Descendants(obj)
            % returns children of children tubes in a breadth-first search
            % list inclusive.
            %
            % .. note:
            %   if tubes can have multiple parents there may be
            %   duplicates.
            %
            descendants = obj;
            tubes = obj;
            while ~isempty(tubes)
                currenttube = tubes(1);
                descendants = [descendants currenttube.children];
                tubes = [tubes currenttube.children];
                tubes(1) = [];
            end
        end

        function ancestors = Ancestors(obj)
            % returns parents of parent tubes inclusive.
            %
            % .. note:
            %   if tubes have multiple parents there may be
            %   duplicates.
            %
            ancestors = obj;
            tubes = obj;
            while ~isempty(tubes)
                currenttube = tubes(1);
                ancestors = [ancestors currenttube.parent];
                tubes = [tubes currenttube.parent];
                tubes(1) = [];
            end
        end

        % spline related
        function obj = MakeSpline(obj, options)
            % fits a spline to the centreline of the tube.
            %
            % Using the object property :attr:`tube.Tube.skelpoints`, a polynomial
            % spline is fit to this list of points. A moving average is
            % used to smooth the skeletal points, it can also use the
            % parent tube to initialise this moving average. For a better
            % understanding of the spline output see `cscvn`_.
            %
            % .. note:
            %   the order of skeletal points affects outcome e.g. reversing
            %   the order of the skeleton points would reverse the spline
            %   direction.
            %
            % Args:
            %   useparent(bool): *OPTIONAL* `default = true` uses the parent tube skel
            %       points if available to initialise the moving average.
            %
            %
            arguments
                obj
                options.useparent logical = true
            end

            % get linear indexed points of previous branch if available.
            if options.useparent == true && isfield(obj.relatives,'parent')
                parent_points = obj.relatives.parent.skelpoints;
                [x_p1, y_p1, z_p1] = I2S(size(obj.network.source), parent_points);
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

        function obj = FindSplinePoints(obj, options)
            % find equidistant points at some given interval along the length of
            % the spline.
            %
            % Integrates along the discontinuous portions of the spline at the
            % chosen interval to save equidistant parametrized points. This is
            % necessary for a continous sampling along the spline otherwise
            % there is no spatial understanding of how often we are resampling
            % a tube. e.g. curved parts could end up being sampled more.
            % Also calls :meth:`tube.Tube.ComputeTortuosity`.
            %
            %
            % Args:
            %   options.sample_interval(float): *OPTIONAL*
            %       `default = `:attr:`network.TubeNetwork.spline_sample_sz`.
            %        The interval size along which to interpolate the spline.
            %
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
            % Interpolates perpendicular slices as if travelling along the tube
            % spline.
            %
            % Interpolates perpendicular slices along the spline at :attr:`tube.Tube.patchprop`.arcpoints.
            % i.e. the spline interval as found by :meth:`tube.Tube.FindSplinePoints` and saves them to :attr:`tube.Tube.source`.
            %
            % .. note::
            %   This function will use the GPU if available only when
            %   :param:`options.method` = 'linear' and the `parallel copmputing toolbox`_
            %   is installed.
            %
            % .. todo::
            %   * Heavily reliant on the network class structure.
            %       Consider decoupling network in this function.
            %   * Consider incorporating matlabs obliqueslice function.
            %
            % Args:
            %   vol(float): must be three dimensional image
            %   options.type(char): *OPTIONAL* `default = 'infer'`, must be either {'source','seg','infer'}. What type is :param:`vol`.
            %   options.usesegcrop(logical): *OPTIONAL* `default = false` use dynamic slicing, vary the size of the plane dependant on
            %       :attr:`tube.Tube.network.seg` size at that point.
            %   options.method(char): *OPTIONAL* `default = 'cubic'` see `interp3`_ for details.
            %   options.sample_sz = *OPTIOPNAL* `default =
            %       obj.network.plane_sample_sz`. interpolation pixel size. 
            %
            %
            % .. _parallel copmputing toolbox: https://www.mathworks.com/products/parallel-computing.html
            % .. _interp3: https://www.mathworks.com/help/matlab/ref/interp3.html
            %
            arguments
                obj
                vol (:,:,:)
                options.type char {mustBeMember(options.type,{'source','seg','infer'})} = 'infer'
                options.usesegcrop logical = false
                options.method char = 'cubic'
                options.sample_sz = obj.network.plane_sample_sz
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
                    offgrid_val=0, method=options.method);
            end

            % Save plane images
            if strcmp(reformedproperty, 'source')
                obj.source = reformedimages;
            elseif strcmp(reformedproperty, 'seg')
                obj.seg = reformedimages;
            end
        end

        function segdiameter = ApproxSegDiameter(obj, sourcepoint)
            % Gets the approximate diameter of the segmentation at a given point.
            %
            % Uses the distance transform of :attr:`tube.Tube.network.seg`
            %   to get the approximate diameter of the segmentation at a given point.
            %
            % .. todo:: * add documentation to this function
            %   * consider diameter conversion to mm
            %
            % Args:
            %   sourcepoint(int): triplet coordinate of the point to sample.
            %
            % Return:
            %   1 variable
            %   * segdiameter(float): segmentation diameter at that point.
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
            % Computes euclidean length of the tube.
            %
            % Computes euclidean length of the tube by considering the skeleton
            % endpoints and saves to :attr:`tube.Tube.stats`.euclength
            %
            %

            assert(isfield(obj.patchprop,'parapoints'), 'no parapoints computed, see method ComputeSplinePoints')
            [~, point_1] = spline_normal(obj.spline, ...
                obj.patchprop.parapoints(1));
            [~, point_end] = spline_normal(obj.spline, ...
                obj.patchprop.parapoints(end));
            obj.stats.euclength = norm(point_end - point_1);
        end

        function obj = ComputeArcLength(obj)
            % Computes arclength of the tube.
            %
            % Computes arc length of the tube by considering the total length
            % of the tube spline and saves to :attr:`tube.Tube.stats`.arclength.
            %
            %

            obj.stats.arclength = Compute_Spline_Points(obj.spline);
        end

        function obj = ComputeTortuosity(obj)
            % Computes tortuosity of the tube.
            %
            % Computes tortuosity of the tube by considering the ratio of the
            % arc:euclidean length and saves to :attr:`tube.Tube.stats`.tortuosity.
            %
            %
            %

            if ~isfield(obj.stats,'arclength')
                obj = ComputeArcLength(obj);
            end
            obj = ComputeEucLength(obj);
            % arclength / euclidean length
            obj.stats.tortuosity = obj.stats.arclength./obj.stats.euclength;
            assert(obj.stats.tortuosity >= 1, 'Impossible to get a tortuosity > 1')
        end

        function meanDval = ComputeMeanDiameter(obj, trim)
            % Computes the trim mean diameter of each measurement type in
            % :attr:`tube.Tube.diameters`.
            %
            % Where :attr:`tube.Tube.diameters` is a `m x n` matrix, the trimmean
            % is calculated for each m row returning an n length vector. result is saved
            % in `tube.Tube.stats`.trimmean and  `tube.Tube.stats`.trimmean_trim.
            %
            %
            % Args:
            %   trim(float): *OPTIONAL* `default = 0` trim as % of extremes of d
            %       ata to discard in trimmean calculation.
            %
            %

            if nargin < 2
                trim = 0;
            end

            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')

            % prune the two variables
            pruned = obj.PruneMeasure(obj.diameters);

            % compute average
            meanDval = trimmean(pruned', trim);
            obj.stats.trimmean = meanDval;
            obj.stats.trimmean_trim = trim;
        end

        function meanAval = ComputeMeanArea(obj, trim)
            % Computes the mean area of each measurement type in
            % :attr:`tube.Tube.areas`.
            %
            % Where :attr:`tube.Tube.areas` is a `m x n` matrix, the mean
            % is calculated for each m row returning an n length vector. result is saved
            % in `tube.Tube.stats`.meanarea and `tube.Tube.stats`.meanarea_trim`.
            %
            % Args:
            %   trim(float): *OPTIONAL* `default = 0` trim as % of extremes of
            %       data to discard in trimmean calculation.
            %
            if nargin < 2
                trim = 0;
            end

            assert(~isempty(obj.areas), 'No areas property. Need measurements.')

            % prune the two variables
            pruned = obj.PruneMeasure(obj.areas);

            % compute average
            meanAval = trimmean(pruned', trim);
            obj.stats.meanarea = meanAval;
            obj.stats.meanarea_trim = trim;
        end

        function intrataperval = ComputeIntrataper(obj)
            % Computes the intrataper value of the tube.
            %
            % the intrataper is the % of tapering in diameter along the length
            % of the tube. It is calculated by using a bisquare linear fitting
            % and then taking the ratio of the coefficients of the result.
            %
            % The result is saved in `tube.Tube.stats`/intrataper.
            %
            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')
            % prune the two variables
            al = obj.PruneMeasure(obj.patchprop.arcpoints);
            var = obj.PruneMeasure(obj.diameters);
            % fit bisquare method
            nrings = size(obj.diameters,1);
            coeff = NaN(nrings,2);
            for ii = 1:nrings
                coeff(ii,:) = robustfit(al, var(ii,:),'bisquare');
            end
            % compute intra-branch tapering as percentage
            intrataperval = -coeff(:,2)./coeff(:,1) * 100;
            obj.stats.intrataper = intrataperval;
        end

        function gradientval = ComputeGradient(obj)
            % Computes the gradient value of the tube.
            %
            % the gradient is the % gradient of tapering in diameter along the length
            % of the tube. It is calculated by using a bisquare linear fitting
            % and then taking gradient coefficient of the result.
            %
            % The result is saved in `tube.Tube.stats`.gradient.
            %
            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')
            % prune the two variables
            al = obj.PruneMeasure(obj.patchprop.arcpoints);
            var = obj.PruneMeasure(obj.diameters);
            % fit bisquare method
            nrings = size(obj.diameters,1);
            coeff = NaN(nrings,2);
            for ii = 1:nrings
                coeff(ii,:) = robustfit(al, var(ii,:),'bisquare');
            end
            % compute intra-branch tapering as percentage
            gradientval = -coeff(:,2) * 100;
            obj.stats.gradient = gradientval;
        end

        function intertaperval = ComputeIntertaper(obj, trim)
            % Computes the intertaper value of the tube.
            %
            % the intertaper value is the % change in mean diameter
            % of the tube relative to the parent tube's mean diameter.
            %
            % The result is saved in `tube.Tube.stats`.intertaper.
            %
            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')
            assert(length(obj.parent) < 2, ['Must only have max one ' ...
                'parent tube. Got ',num2str(length(obj.parent))])

            % compute means
            parentmean = obj.parent.ComputeMeanDiameter(trim);
            currentmean = obj.ComputeMeanDiameter(trim);

            % compute interbranch tapering as percentage
            intertaperval = (parentmean - currentmean)./(parentmean) * 100;
            obj.stats.intertaper = intertaperval;
        end

        function volumeval = ComputeVolume(obj)
            % Computes the volume of the tube.
            %
            % the volume is calculated as an integration of the diameters against
            % its arclength.
            %
            % The result is saved in `tube.Tube.stats`.volume.
            %
            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')
            % prune the two variables
            al = obj.PruneMeasure(obj.patchprop.arcpoints);
            var = obj.PruneMeasure(obj.diameters);
            % fit bisquare method
            nrings = size(obj.diameters,1);
            volumeval = NaN(nrings,1);
            for ii = 1:nrings
                volumeval(ii) = trapz(al, var(ii,:));
            end
            % compute intra-branch tapering as percentage
            obj.stats.volume = volumeval;
        end

        % measures
        function obj = Measure(obj, classmethod, varargin)
            % Call desired method to make measurement on the interpolated tube slices.
            %
            % This calls a subclass of the superclass, :class:`measure.SuperMeasure`
            % to make measurement on the interpolated slices of the tube.
            %
            % .. todo:: add documentation to this function
            %
            % Args:
            %   classmethod(char):
            %   varargin
            %
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

        % 2D visualisation
        function h = plot(obj, X, Y)
            arguments
            obj
            X = obj.patchprop.arcpoints
            Y = obj.diameters
            end
            % plot patchprop measure
            %
            % desc
            %
            % .. todo::
            %   * add documentation to this function
            %   * Needs attention
            %
            % Args:
            %  relativetube (:class:`tube`): the tube to set
            %   relation to.
            %  relation (string): relation name. common
            %   "parent" or "child".
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.tubes(98).plot();
            %
            % .. |tube_Tube_plot| image:: figs/tube_plot.png
            %    :width: 400
            %    :alt: figure plot - Tube plot
            %

            X = parsearg(X);
            Y = parsearg(Y);

            assert(all(size(X,2) == size(Y,2)), ['X and Y needs have same n cols' ...
                '. Got X [', num2str(size(X,2)), ['] and Y [' ...
                ''], num2str(size(Y,2)),']'])

            X = repmat(X, size(Y,1), 1);

            h = plot(X', Y', '.');
            legend(arrayfun(@(mode) sprintf('Ring %d', mode), 1:size(Y, 2), 'UniformOutput', false))

            % default color order
            colororder(linspecer(12))

            function parsedarg = parsearg(arg)
                if isa(arg, "char") || isa(arg, "string")
                    parsedarg = obj.(arg);
                else
                    parsedarg = arg;
                end
            end
        end

        function s = OrthoView(obj, options)
            % View a series of an airway segment's slices as a volume image
            % stack using MATLAB's inbuilt othogonal 3d viewer.
            % short desc
            %
            % long desc
            %
            % .. todo::
            %   * add documentation to this function
            %   * add 4th panel to show plot
            %   * make scrollable by mousewheel
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.tubes(98).OrthoView();
            %
            % .. |tube_Tube_OrthoView| image:: figs/tube_orthoview.png
            %    :width: 400
            %    :alt: figure plot - Tube Orthoview
            %
            % |tube_Tube_OrthoView|
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
            tubearray = ParseVolOut(obj, type=options.type);

            % display with orthoview
            s = orthosliceViewer(tubearray, 'DisplayRangeInteraction','off', ...
                'ScaleFactors',obj.voxdim, 'CrosshairLineWidth', 0.3);
            % get axial axes
            [ax, ~, ~] = getAxesHandles(s);

            addlistener(s,'CrosshairMoving',@allevents);
            addlistener(s,'CrosshairMoved',@allevents);

            % default color order
            colororder(linspecer(12))

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

        % 3D visualisation - level 1
        function h = Plot3(obj, color)
            % plot this tube as a line asif an edge on a 3D network graph.
            %
            % .. todo: consider adding more info in datatips.
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.tubes(98).Plot3();
            %
            % .. |tube_Tube_Plot3| image:: figs/tube_plot3.png
            %    :width: 400
            %    :alt: figure plot - Tube Plot3
            %
            % |tube_Tube_Plot3|
            %
            if nargin < 2
                color =[];
            end
            % get endpoints
            [Y,X,Z] = obj.I2S([obj.skelpoints(1),obj.skelpoints(end)]);
            h = plot3(X,Y,Z,'.');
            h.Color = 'k';
            hold on
            % get midpoint point of line
            [mx,my,mz] = midpoint(X,Y,Z);
            X = [X(1), mx, X(2)];
            Y = [Y(1), my, Y(2)];
            Z = [Z(1), mz, Z(2)];
            % plot line
            h = plot3(X,Y,Z,'-');
            % datatip

            dtRows = [dataTipTextRow("tubeID",ones(1,3)*obj.ID), ...
                dataTipTextRow("generation",ones(1,3)*obj.generation)];

            h.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows;

            % colour code?
            if ~isempty(color)
                h.Color = color;
            end
        end

        function h = Plot3D(obj, options)
            %
            %
            % .. todo: add datatip for stats
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.tubes(98).Plot3D();
            %
            % .. |tube_Tube_Plot3D| image:: figs/tube_plot3d.png
            %    :width: 400
            %    :alt: figure plot - Tube Plot3D
            %
            % |tube_Tube_Plot3D|
            %

            arguments
                obj
                options.type {mustBeMember(options.type,{'seg','skel'})} = 'seg'
                options.alpha = 0.3
                options.color = 'k'
                options.context = true;
                options.contextcolor = 'r'
            end
            if options.context == true
                % get adjacent tubes
                adjtubes = [obj.parent, obj.children, obj.parent(1).children];
                for adjii = adjtubes
                    adjii.Plot3D(type=options.type, color=options.contextcolor,...
                        context=false, alpha=options.alpha)
                    hold on
                end
            end
            V = zeros(size(obj.network.seg));

            V(obj.([options.type, 'points'])) = 1;

            h = patch(isosurface(V),...
                'FaceAlpha', options.alpha,...
                'FaceColor', options.color,...
                'EdgeColor', 'none');
            obj.network.vol3daxes()


            hold off
        end

        function PlotSpline(obj,options)
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.tubes(98).PlotSpline();
            %
            % .. |tube_Tube_PlotSpline| image:: figs/tube_plotspline.png
            %    :width: 400
            %    :alt: figure plot - Tube PlotSpline
            %
            % |tube_Tube_PlotSpline|
            %
            arguments
                obj
                options.color = 'k'
                options.context = true;
                options.contextcolor = 'r'
            end

            if options.context == true
                % get adjacent tubes
                adjtubes = [obj.parent, obj.children, obj.parent(1).children];
                for adjii = adjtubes
                    adjii.PlotSpline(color=options.contextcolor,...
                        context=false)
                    hold on
                end
            end

            rawpoints = fnplt(obj.spline);
            points = rawpoints./obj.network.voxdim';
            h = plot3(points(2,:),points(1,:),points(3,:),'Color',options.color);

            dtRows = [dataTipTextRow("tubeID",ones(size(points,2))*obj.ID), ...
                dataTipTextRow("generation",ones(size(points,2))*obj.generation)];

            h.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows;

            obj.network.vol3daxes();
            hold off
        end

        % 3D visualisation - level 2
        function h = PlotSplineVecs(obj, options)
            %
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.tubes(98).PlotSplineVecs();
            %
            % .. |tube_Tube_PlotSplineVecs| image:: figs/tube_plotsplinevecs.png
            %    :width: 400
            %    :alt: figure plot - Tube PlotSplineVecs
            %
            % |tube_Tube_PlotSplineVecs|
            %
            arguments
            obj
            options.subsamp = 2
            end

            % viewvol for adjacent tubes
            obj.Plot3D()
            hold on

            % get spline points and their normals
            % subsample data, for better visualisation
            samplepnts = obj.patchprop.parapoints(1:options.subsamp:end);
            origins = NaN(3,size(samplepnts,2));
            normvecs = NaN(3,size(samplepnts,2));

            for ii = 1:length(samplepnts)
                [normvecs(:,ii), origins(:,ii)] = spline_normal(obj.spline, ...
                    samplepnts(ii));
            end

            % rescale origins from mm to vox and swap x and y in plot.
            origins = origins./obj.network.voxdim';
            normvecs = normvecs./obj.network.voxdim';

            % plot vectors
            h = quiver3(origins(2,:),origins(1,:),origins(3,:),...
                normvecs(2,:),normvecs(1,:),normvecs(3,:),'k');
            hold off

            obj.network.vol3daxes()
        end

        % Data IO

        function obj = ExportPerpPathces(obj,path)
            % export perpendicular slice patches of this tube.
            %
            % export the perpendicular slice patches of this tube stored in
            % source as int16 tiff files.
            %
            %
            % Args:
            %   path(str): path to directory to save the exported patches.
            %       The directory will be created if it doesn't already 
            %       exist.
            %
            %

            % make directory
            dirname = fullfile(path,'airway_patches');
            if ~exist(dirname, 'dir')
                mkdir(dirname)
            end

            % choose which slices to save
            chosenslices = PruneMeasure(obj, 1:length(obj.source));
            % loop through slices
            for k = chosenslices
                img = int16(obj.source(:,:,k));

                % save as int16 TIF
                imgsavename = fullfile(dirname, [ ...
                    saveid, '_', ...
                    'id_',num2str(obj.id), ...
                    '_gen_', num2str(obj.generation), ...
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
            end
        end

        function toGif(obj, filename, options)
            % View a series of an airway segment's slices as a volume image
            % stack using MATLAB's inbuilt othogonal 3d viewer.
            % short desc
            %
            % long desc
            %
            % .. todo:: add documentation to this function
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
                options.framerate (1,1) mustBeNumeric = 20
            end

            % parse filename
            filename = parse_filename_extension(filename, '.gif');

            % instantiate orthosliceviewer
            s = obj.OrthoView(type=options.type, rings=options.rings, ...
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
                    imwrite(indI,cm,filename,'gif','Loopcount',inf,...
                        'DelayTime',1/options.framerate);
                else
                    imwrite(indI,cm,filename,'gif','WriteMode',...
                        'append','DelayTime',1/options.framerate);
                end
            end
        end

        % utilities
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

            I = S2I3(size(obj.network.source),I1,I2,I3);
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
            % var is of dimensions nrings x nslices
            if isempty(obj.prunelength)
                obj.prunelength = [0 0];
                warning('Object prunelength not set, using zero pruning.')
            end
            nrings = size(var,1);
            pl = obj.prunelength;
            al = obj.patchprop.arcpoints;
            assert(length(al) == length(var), 'Input variable must be the same length as arclength.')

            prunebool = (al >= pl(1) & al <= (al(end) - pl(2)));
            % scale to number of measures
            prunebool = repmat(prunebool,nrings,1);
            prunedprop = var(prunebool == 1);
            prunedprop = reshape(prunedprop,nrings,[]);
        end

        function volout = ParseVolOut(obj,options)
            % short desc
            %
            % long desc
            %
            % .. todo::
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
            volout = zeros([canvas_sz, canvas_sz, length(tubecell)]);
            for slice = 1:length(tubecell)
                image = tubecell{slice,1};
                image_sz = size(image,1);
                min_centre = canvas_sz/2 - image_sz/2;
                max_centre = canvas_sz/2 + image_sz/2;
                volout(min_centre+1:max_centre, min_centre+1:max_centre, slice) = image;
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
