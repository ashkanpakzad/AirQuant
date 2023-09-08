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
        sibling = []
        generation
        network
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
        lead_dir = []
        trail_dir = []
    end

    methods
        % Core
        function obj = Tube(network, skelpoints, ID)
            % Init Tube see :class:`tube.Tube`.
            %
            % Makes the Tube object by defining the skelpoints along the centreline of the tube.
            % This init method automatically interpolates the tube and fits a spline.
            % The network object is required as it specifies certain configuration properties.
            % This function is not intended to be invoked directrly. Only by
            % :meth:`network.TubeNetwork.MakeTubes`.
            %
            % .. todo::
            %   * Generalise this init function so that empty tubes can be
            %       created.
            %   * make it possible for tubes to be created independent of
            %       network object.
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

            obj.patchprop = struct;
            obj.stats = struct;
            obj.region = struct;

            obj.MakeSpline();
            obj.FindSplinePoints();

        end
        
        function obj = SetRoot(obj)
            % Set this tube as the root in the network.
            %
            % The root tube is defined to have no parents and that every connected tube 
            % will be a decendent. This will manipulate every connected 
            % tube to ensure it originates from this tube. Tubes which are
            % flipped in relation (i.e. parent becomes child and vice
            % versa) will also have its skel points flipped and spline 
            % remade.
            %

            if isempty(obj.parent)
                % already root
                return
            end
            root = obj;

            % init variables
            newparent = obj;
            oldparent = newparent.parent;
            % loop through all parents
            while ~isempty(oldparent)
                % old siblings become children
                siblings = newparent.GetSibling();
                for sibling = siblings
                    % eliminate direct relations with old parent
                    sibling.parent(find(sibling.parent == oldparent)) = [];
                    oldparent.children(find(oldparent.children == sibling)) = [];
                    % set sibling as child
                    newparent.SetChildren(sibling);
                end
                
                % flip parent and children roles
                newparent.SetChildren(oldparent)
                % delete old references
                oldparent.children(find(oldparent.children == newparent)) = [];
                newparent.parent(find(newparent.parent == oldparent)) = [];
                
                % flip parent indices
                newparent.skelpoints = flip(newparent.skelpoints);
                % redo spline with flipped indices
                newparent.MakeSpline();
                newparent.FindSplinePoints();
                
                % prepare next iteration
                newgrandchild = oldparent.parent(find(oldparent.parent ~= newparent));
                newparent = oldparent;
                oldparent = newgrandchild;
            end

            % set ID by BFS
            root.SetID_BFS()
        end

        function obj = SetID_BFS(obj)
            % set ID of all tubes by BFS from current tube.
            %
            % .. warning:
            %   This will reset all the IDs of tubes. IDs will be duplicates
            %   if the tube network is not a complete component.
            %
            id = 1;
            obj.ID = id;
            todo = obj.children;
            while ~isempty(todo)
                current = todo(1);
                todo(1) = [];
                id = id + 1;
                current.ID = id;
                todo = [todo, current.children];
            end

            % sort array by new id in tube network object
            [~,idx]=sort([obj.network.tubes.ID]);
            obj.network.tubes = obj.network.tubes(idx);
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

        function siblings = GetSibling(obj)
            % Get siblings of tube if they exist.
            %
            %
            % Returns:
            %   1 variable.
            %   * siblings(`vector of tube`) = returns vector of tubes. If
            %     no parent or no siblings then returns empty.
            %
            
            % check if parent exists
            if isempty(obj.parent)
                siblings = [];
                return
            end
            
            parent_children = obj.parent.children;
            % check if any siblings
            if length(parent_children) == 1
                siblings = [];
                return
            end
            
            % remove current tube to get siblings
            siblings = parent_children(parent_children ~= obj);
        end

        function obj = SetGeneration(obj)
            % Identify and set the generation of this tube.
            %
            % Counts the number of tube anscestors i.e. parents of parents up to
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
        % check same region
%         if length(currentube.parent) > 1
        while ~isempty(currentube.parent) 
            current_reg = currentube.region.(regiontype);
            parent_reg = currentube.parent.region.(regiontype);
            if ~strcmp(current_reg, parent_reg)
                break
            end
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
            if options.useparent == true && ~isempty(obj.parent)
                parent_points = obj.parent.skelpoints;
                [x_p1, y_p1, z_p1] = I2S(obj,parent_points);
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
            [arclength, obj.patchprop.parapoints, ...
                obj.patchprop.arcpoints] = ...
                Compute_Spline_Points(obj.spline, options.sample_interval);

            obj.stats.arclength = real(arclength);
            % save stats measurement using derived spline points.
            obj = ComputeTortuosity(obj);
        end
        
        function [lead_dir, trail_dir] = ComputeDirections(obj)
            % Compute the leading and trailing directions of the tube.
            %
            % Using the beginning and end section of the tube, calculates
            % the leading and trailing direction vectors respectively.
            %
            % Returns:
            %   2 variables
            %   * lead_dir(`3x1 vector`) = vector of leading direction.
            %   * trail_dir(`3x1 vector`) = vector of trailing direction.
            %

            % compute leading direction from beginning spline
            % where spline parameter t = 0
            lead_dir = spline_normal(obj.spline, 0);
            obj.lead_dir = lead_dir;
            % compute trailing direction from end of spline
            % where spline parameter t = end
            trail_dir = spline_normal(obj.spline, obj.spline.breaks(end));
            obj.trail_dir = trail_dir;

        end

        function [change_deg] = ComputeChangeAngle(obj)
            % Compute the change in angle this tube.
            %
            % Compute the change in between the leading direction to the
            % trailing direction of this tube.
            %
            % Returns:
            %   1 variables
            %   * change_deg(`scalar`) = angle in degrees. 
            %

            P1 = obj.lead_dir;
            P2 = obj.trail_dir;
            change_deg = vector_angle(P1,P2,'degrees');

            obj.stats.change_deg = change_deg;
        end

        function [parent_deg] = ComputeParentAngle(obj)
            % Compute the angle of the parent to this tube.
            %
            % Compute the angle of the parent trailing direction to the
            % leading direction of this tube. If this tube has no parent
            % then the output will be NaN.
            %
            % Returns:
            %   1 variables
            %   * parent_deg(`scalar`) = angle in degrees. NaN if no parent.
            %

            % check if parent exists
            if isempty(obj.parent)
                parent_deg = NaN;
            else
            P1 = obj.parent.trail_dir;
            P2 = obj.lead_dir;
            parent_deg = vector_angle(P1,P2,'degrees');
            end

            obj.stats.parent_deg = parent_deg;
        end

        function [sibling_deg] = ComputeSiblingAngle(obj,overwrite)
            % Compute the angle of this tube to its sibling (bifurcation
            % only).
            %
            % Computes the leading angle of this tube to the
            % leading direction of its sibling. If this tube has multiple
            % siblings then the output will be NaN.
            %
            % Args:
            %   overwrite(`bool`) = OPTIONAL `default = true`. If false and
            %   value has already been set then completes. This behaviour
            %   has been added for efficiency.
            %
            % Returns:
            %   1 variables
            %   * sibling_deg(`scalar`) = angle in degrees. NaN if multiple
            %       siblings.
            %
            
            if nargin < 2
                overwrite = 1;
            end

            % only complete if not already set and overwrite is off
            if isfield(obj.stats,'sibling_deg') && overwrite == 0
                sibling_deg = [];
                return
            end

            % must be only 1 sibling.
            siblings = obj.GetSibling();
            if length(siblings) ~= 1
                sibling_deg = NaN;
            else
                P1 = siblings.lead_dir;
                P2 = obj.lead_dir;
                sibling_deg = vector_angle(P1,P2,'degrees');
                % save in sibling tube
                siblings.stats.sibling_deg = sibling_deg;
            end
            obj.stats.sibling_deg = sibling_deg;
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
            %   * GPU only compatible when :param:`options.method` = 'linear' and the `parallel copmputing toolbox`_
            %       is installed.
            %   * GPU is not compatible with :param:`options.usesegcrop` = true.
            %   * GPU is a legacy feature and will be removed in future versions.
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
                options.method char = 'linear'
                options.sample_sz = obj.network.plane_sample_sz
                options.gpu logical = false
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
                    sample_sz=options.sample_sz, ...
                    offgrid_val=0, method=options.method, gpu=options.gpu);
            end

            % Save plane images
            if strcmp(reformedproperty, 'source')
                obj.source = reformedimages;
            elseif strcmp(reformedproperty, 'seg')
                obj.seg = reformedimages;
            end
        end

        function segdiameter = ApproxSegDiameter(obj, sourcepoint)
            % Gets the approximate diameter using the segmentation for a given point.
            %
            % Uses the distance transform of :attr:`tube.Tube.network.seg`
            %   to get the approximate diameter of the segmentation at a given point.
            %
            % .. todo:: 
            %   * consider diameter conversion to mm
            %
            % Args:
            %   sourcepoint(int,triplet): subindex coordinates of the point to sample.
            %
            % Return:
            %   1 variable
            %   * segdiameter(float): segmentation diameter in pixels.
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
        function euclength = ComputeEucLength(obj)
            % Computes euclidean length of the tube.
            %
            % Computes euclidean length of the tube by considering the skeleton
            % endpoints and saves to :attr:`tube.Tube.stats`.euclength
            %
            % .. note:
            %   Requires splines to have been fitted first using
            %   :meth:`tube.Tube.FindSplinePoints`.
            %

            assert(isfield(obj.patchprop,'parapoints'), 'no parapoints computed, see method ComputeSplinePoints')
            [~, point_1] = spline_normal(obj.spline, ...
                obj.patchprop.parapoints(1));
            [~, point_end] = spline_normal(obj.spline, ...
                obj.patchprop.parapoints(end));
            euclength = real(norm(point_end - point_1));
            if euclength == 0
                warning('Euclidean length of a branch is 0, try reducing your spline sampling size')
            end
            obj.stats.euclength = euclength;
        end

        function obj = ComputeArcLength(obj)
            % Computes arclength of the tube.
            %
            % Computes arc length of the tube by considering the total length
            % of the tube spline and saves to :attr:`tube.Tube.stats`.arclength.
            %
            % .. note:
            %   * Requires splines to have been fitted first using
            %   :meth:`tube.Tube.FindSplinePoints`.
            %   * This function is not called by class internals,
            %   :meth:`tube.Tube.FindSplinePoints` does instead to avoid
            %   code redundancy.

            obj.stats.arclength = real(Compute_Spline_Points(obj.spline));
        end

        function tortuosity = ComputeTortuosity(obj)
            % Computes tortuosity of the tube.
            %
            % Computes tortuosity of the tube by considering the ratio of the
            % arc:euclidean length and saves to :attr:`tube.Tube.stats`.tortuosity.
            %
            % .. note:
            %   Requires splines to have been fitted first using
            %   :meth:`tube.Tube.ComputeSplinePoints.
            %

            if ~isfield(obj.stats,'arclength')
                obj = ComputeArcLength(obj);
            end
            ComputeEucLength(obj);
            % arclength / euclidean length
            tortuosity = real(obj.stats.arclength./obj.stats.euclength);
            % note that we evaluate with a tolerance for the rare case that
            % arclength is very close to 1 due to the precision of
            % numerical methods.
            precision = 1e-6;
            if tortuosity < 1-precision
                warning(strcat("Tortuosity shouldn't be less than 1 " + ...
                    "by definition. Got ", ...
                num2str(obj.stats.tortuosity),". This could be " + ...
                    "due to numerical imprecision."))
            end

            obj.stats.tortuosity = tortuosity;
        end

        function [meanDminval, meanDmajval] = ComputeMeanMinMajDiameters(obj, trim)
            % Computes the trim mean diameter of the minor and major diameters of each measurement type in
            % :attr:`tube.Tube.diameters`.
            %
            % Where :attr:`tube.Tube.diameters` is a `m x n` matrix, the trimmean
            % is calculated for each m row returning an n length vector. result is saved
            % in `tube.Tube.stats`.trimmean and  `tube.Tube.stats`.trimmean_trim. 
            % For ordinary mean (no trim) set trim to 0.
            %
            % Args:
            %   trim(float): *OPTIONAL* `default = 0` trim as % of extremes 
            %   of data to discard in trimmean calculation.
            %
            %

            if nargin < 2
                trim = 0;
            end

            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')        

            % prune the two variables
            pruned_min = obj.PruneMeasure(obj.patchprop.min_diameters);
            pruned_maj = obj.PruneMeasure(obj.patchprop.maj_diameters);

            % compute average
            meanDminval = real(trimmean(pruned_min', trim));
            meanDmajval = real(trimmean(pruned_maj', trim));
            obj.stats.min_diameter_mean = meanDminval;
            obj.stats.maj_diameter_mean = meanDmajval;
            obj.stats.minmajdiameter_trim = trim;
        end

        function meanDval = ComputeMeanDiameter(obj, trim)
            % Computes the trim mean diameter of each measurement type in
            % :attr:`tube.Tube.diameters`.
            %
            % Where :attr:`tube.Tube.diameters` is a `m x n` matrix, the trimmean
            % is calculated for each m row returning an n length vector. result is saved
            % in `tube.Tube.stats`.trimmean and  `tube.Tube.stats`.trimmean_trim. 
            % For ordinary mean (no trim) set trim to 0.
            %
            % Args:
            %   trim(float): *OPTIONAL* `default = 0` trim as % of extremes 
            %   of data to discard in trimmean calculation.
            %
            %

            if nargin < 2
                trim = 0;
            end

            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')

            % prune the two variables
            pruned = obj.PruneMeasure(obj.diameters);

            % compute average
            meanDval = real(trimmean(pruned', trim));
            obj.stats.diameter_mean = meanDval;
            obj.stats.diameter_trim = trim;
        end

        function meanPerimeter= ComputeMeanPerimeter(obj, trim)
            % Computes the trim mean perimeter of each measurement type in
            % :attr:`tube.Tube.patchprop.perimeter.
            %
            % Args:
            %   trim(float): *OPTIONAL* `default = 0` trim as % of extremes
            %   of data to discard in trimmean calculation.
            %
            %

            if nargin < 2
                trim = 0;
            end

            assert(~isempty(obj.patchprop.perimeter), ...
                'No perimeter in patchprop. Need measurements.')

            % prune the two variables
            pruned = obj.PruneMeasure(obj.patchprop.perimeter);

            % compute average
            meanPerimeter = real(trimmean(pruned', trim));
            obj.stats.perimeter_mean = meanPerimeter;
            obj.stats.perimeter_trim = trim;
        end


        function meanHydraulicDval = ComputeMeanHydraulicD(obj, trim)
            % Computes the trim mean hydraulic diameter of each measurement type in
            % :attr:`tube.Tube.patchprop.hydraulic_diameter`.
            %
            % Where :attr:`tube.Tube.patchprop.hydraulic_diameter` is a `m x n` matrix, the trimmean
            % is calculated for each m row returning an n length vector. result is saved
            % in `tube.Tube.stats`.hydraulicD_mean and  `tube.Tube.stats`.hydraulicD_trim.
            % For ordinary mean (no trim) set trim to 0.
            %
            % Args:
            %   trim(float): *OPTIONAL* `default = 0` trim as % of extremes
            %   of data to discard in trimmean calculation.
            %
            %

            if nargin < 2
                trim = 0;
            end

            assert(~isempty(obj.patchprop.hydraulic_diameter), ...
                'No hydraulic_diameter in patchprop. Need measurements.')

            % prune the two variables
            pruned = obj.PruneMeasure(obj.patchprop.hydraulic_diameter);

            % compute average
            meanHydraulicDval = real(trimmean(pruned', trim));
            obj.stats.hydraulicD_mean = meanHydraulicDval;
            obj.stats.hydraulicD_trim = trim;
        end

        function meanAval = ComputeMeanArea(obj, trim)
            % Computes the mean area of each measurement type in
            % :attr:`tube.Tube.areas`.
            %
            % Where :attr:`tube.Tube.areas` is a `m x n` matrix, the mean
            % is calculated for each m row returning an n length vector. result is saved
            % in `tube.Tube.stats`.meanarea and `tube.Tube.stats`.meanarea_trim.
            %
            % Args:
            %   trim(float): *OPTIONAL* `default = 0` trim as % of extremes 
            %       of data to discard in trimmean calculation.
            %
            % Returns:
            %   1 variable
            %   * meanAval (array) mean value of each ring.
            %

            if nargin < 2
                trim = 0;
            end

            assert(~isempty(obj.areas), 'No areas property. Need measurements.')

            % prune the two variables
            pruned = obj.PruneMeasure(obj.areas);

            % compute average
            meanAval = real(trimmean(pruned', trim));
            obj.stats.area_mean = meanAval;
            obj.stats.area_trim = trim;
        end

        function intrataperval = ComputeIntrataper(obj)
            % Computes the intrataper value.
            %
            % Intrataper is the % of tapering in diameter along the length
            % of the tube. It is calculated by using a bisquare linear
            % fitting to fit the line :math:`y=mx+c` to measurements.
            % Then computing the ratio of the coefficients.
            % The result is saved in `tube.Tube.stats`intrataper.
            %
            % .. math::
            %   intrataper = - m / c
            %
            % Returns:
            %   1 variable
            %   * intrataperval (array) intrataper value of each ring.
            %
            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')
            % prune the two variables
            al = obj.PruneMeasure(obj.patchprop.arcpoints);
            var = obj.PruneMeasure(obj.diameters);
            % fit bisquare method
            nrings = size(obj.diameters,1);
            coeff = NaN(nrings,2);
            intrataperval = NaN(nrings, 1);
            try % incase no branch left after pruning/too few points
                for ii = 1:nrings
                    coeff(ii,:) = robustfit(al, var(ii,:),'bisquare');
                    intrataperval(ii) = real(-coeff(ii,2)./coeff(ii,1) * 100);
                end
            catch
                % leave as nan
            end
            % compute intra-branch tapering as percentage
            obj.stats.intrataper = intrataperval;
        end

        function gradientval = ComputeGradient(obj)
            % Computes the gradient value of the tube.
            %
            % the gradient is the % gradient of tapering in diameter along the length
            % of the tube. It is calculated by using a bisquare linear
            % fitting to fit the line :math:`y=mx+c` to measurements.
            % and then taking gradient coefficient of the result. 
            % The result is saved in `tube.Tube.stats`.gradient.
            %
            % .. math::
            %   gradient = -m
            %
            % Returns:
            %   1 variable
            %   * gradientval (array) gradient value of each ring.
            %

            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')
            % prune the two variables
            al = obj.PruneMeasure(obj.patchprop.arcpoints);
            var = obj.PruneMeasure(obj.diameters);
            % fit bisquare method
            nrings = size(obj.diameters,1);
            coeff = NaN(nrings,2);
            for ii = 1:nrings
                try
                    coeff(ii,:) = robustfit(al, var(ii,:),'bisquare');
                catch % incase robust fit fails
                    coeff(ii,:) = NaN;
                end
            end
            % compute intra-branch tapering as percentage
            gradientval = real(-coeff(:,2) * 100);
            obj.stats.gradient = gradientval;
        end

        function intertaperval = ComputeIntertaper(obj, trim, parentfail)
            % Computes the intertaper value.
            %
            % intertaper value is the % change in mean diameter
            % of the tube relative to the parent tube's mean diameter.
            % The result is saved in `tube.Tube.stats`.intertaper.
            %
            % .. math::
            %   intertaper = \hat{d}_{p} - d / \hat{d}_{p}
            %
            % Args:
            % trim(float): *OPTIONAL* `default = 0` trim as % of extremes 
            %   of data to discard in trimmean calculation.
            % parentfail: *OPTIONAL* `default = false` require lack of
            %   parent measurement to raise error.
            %
            % Returns:
            %   1 variable
            %   * intertaper (array) intertaper value of each ring.
            %

            if nargin < 2
                trim = 0;
            end
            if nargin < 3
                parentfail = false;
            end
            
            assert(~isempty(obj.diameters), 'No diameters property. Need measurements.')
            
            if parentfail % fail if not one parent
                assert(length(obj.parent) == 1, ['Must only have only one ' ...
                    'parent tube. Got ',num2str(length(obj.parent))])
            end
            
            try
                % compute means
                parentmean = obj.parent.ComputeMeanDiameter(trim);
                currentmean = obj.ComputeMeanDiameter(trim);
    
                % compute interbranch tapering as percentage
                intertaperval = real((parentmean - currentmean)./(parentmean) * 100);
            catch
                intertaperval = NaN;
            end
            obj.stats.intertaper = intertaperval;

        end

        function volumeval = ComputeVolume(obj)
            % Computes volume.
            %
            % The volume is calculated as an integration of the areas against
            % its arclength by trapezium numerical approximation. The result 
            % is saved in `tube.Tube.stats`.volume.
            %
            % Returns:
            %   1 variable
            %   * volumeval (`array`) volume value of each ring.
            %
            assert(~isempty(obj.areas), 'No areas property. Need measurements.')
            % prune the two variables
            al = obj.PruneMeasure(obj.patchprop.arcpoints);
            var = obj.PruneMeasure(obj.areas);
            nrings = size(obj.areas,1);
            volumeval = NaN(nrings,1);
            for ii = 1:nrings
                vals = var(ii,:);
                if length(vals) == 1 % incase only one point in tube
                    volumeval(ii) = real(al * vals);
                else
                    volumeval(ii) = real(trapz(al, vals));
                end
            end
            obj.stats.volume = volumeval;
        end

        % measures
        function obj = Measure(obj, classmethod, varargin)
            % Call desired method to make measurement on the interpolated tube patch slices.
            %
            % This calls a subclass of :class:`measure.SuperMeasure`
            % to make measurement on the interpolated patch slices of the tube. 
            % The measure type called must be specific to the tube type.
            % e.g. airway for airway images. Also recomputes all
            % measurements.
            %
            % .. todo:: 
            %   * add more detailed documentation to this function.
            %   * make measurements specific page to link to.
            %   * initiate by deleting all calculations.
            %
            % Args:
            %   classmethod(char): Name of class method to make
            %       measurement.
            %   varargin
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
            % set raw measures to tube properties
            obj.measures = classmeasures.measures;
            obj.patchprop.min_diameters = classmeasures.OutputVar('min_diameter');
            obj.patchprop.maj_diameters = classmeasures.OutputVar('maj_diameter');  
            obj.diameters = classmeasures.OutputVar('diameter');
            obj.areas = classmeasures.OutputVar('area');
            obj.patchprop.perimeter = classmeasures.OutputVar('perimeter');
            % get hydraulic diameter if exist
            try
                obj.patchprop.hydraulic_diameter = classmeasures.OutputVar('hydraulic_diameter');
            catch
            end
            % derive stats of measurements
            obj.ComputeMeanMinMajDiameters();
            obj.ComputeMeanDiameter();
            obj.ComputeMeanHydraulicD();
            obj.ComputeMeanPerimeter();
            obj.ComputeMeanArea();
            obj.ComputeIntrataper();
            obj.ComputeGradient();
            obj.ComputeVolume();
        end

        % 2D visualisation
        function h = plot(obj, options)
            % Line plot two object properties with linear regression
            % line.
            %
            % A versatile method for plotting two object properties or array with a
            % linear regressed line 'line of best fit'. With smoothing if
            % desired.
            %
            % .. todo::
            %   * add visualisation
            %   * link to plot
            %
            % Args:
            %  X (:class:`tube.Tube` property or array): *OPTIONAL* `default =
            %   :attr:`patchprop`.arcpoints` Horizontal X axis property
            %   name.
            %  Y (:class:`tube.Tube` property or array): *OPTIONAL* `default =
            %   :attr:`diameters`` Vertical Y axis property
            %   name.
            %  smoothing (float): *OPTIONAL* `default = 1e-64` degree of
            %   smoothing of the line plot.
            %  linespec (char): *OPTIONAL* `default = -` linespec to 
            %   specify style. see MATLAB's plot function.
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.tubes(98).plot();
            %
            %
            
            arguments
            obj
            options.X = obj.patchprop.arcpoints
            options.Y = obj.diameters
            options.smoothing = 1e-64
            options.linespec = '-'
            end

            % parse inputs, if string, finds matching property
            X = parsearg(options.X);
            Y = parsearg(options.Y);
        
            assert(all(size(X,2) == size(Y,2)), ['X and Y needs have same n cols' ...
                '. Got X [', num2str(size(X,2)), ['] and Y [' ...
                ''], num2str(size(Y,2)),']'])
            
            % matches x to number of y
            X = repmat(X, size(Y,1), 1);
            
            % smooth Y rows individually
            for ii= 1:size(Y,1)
                Y(ii,:) = smooth(Y(ii,:), options.smoothing);
            end
            
            % plot datapoints
            h = plot(X', Y', options.linespec);
            legend(arrayfun(@(mode) sprintf('Measurement %d', mode), 1:size(Y, 2), 'UniformOutput', false))
            
            % linear bestfit line
            for ii= 1:size(Y,1)
                coefficients = polyfit(X(ii,:), Y(ii,:), 1);
                yFit = polyval(coefficients , X(ii,:));
                % Plot everything.
                hold on; 
                plot(X(ii,:), yFit, 'r-'); % Plot fitted line.
                hold off
            end

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
            % View a series of patches of tube slices as a volume image
            % stack using MATLAB's inbuilt othogonal 3d viewer.
            %
            % .. todo::
            %   * add 4th panel to show plot
            %   * make scrollable by mousewheel
            %   * link to orthosliceviewer
            %
            % Args:
            %   type(char): *OPTIONAL* `default = 'source' object property
            %       that is a 3D array to view.
            %   rings(bool array): *OPTIONAL* `default = 
            %       ones(1,size(:attr:`measures`,1))` same length as number
            %       of measures per patch, i.e. number of rings. For each
            %       ring to visualise. All rings by default.
            %   ellipses(bool): *OPTIONAL* `default = true` to view 
            %       ellipses of measure if available.
            %   points(bool): *OPTIONAL* `default = true` to view 
            %       points of measure if available.
            %
            % Return:
            %   1 variable
            %   * s(`orthosliceViewer`): handle to created graphics
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
            obj.UpdateAnnotateOrthoviewer(tubearray,ax,s.SliceNumbers(3),...
                options.rings,options.ellipses...
                ,options.points)

            function allevents(src,evt)
                evname = evt.EventName;
                switch(evname)
                    case{'CrosshairMoving','CrosshairMoved'}
                        ppos = evt.PreviousPosition(3);
                        cpos = evt.CurrentPosition(3);
                        if ppos ~= cpos
                            obj.UpdateAnnotateOrthoviewer(tubearray,ax,cpos,...
                                options.rings,options.ellipses...
                                ,options.points)
                        end
                end
            end

        end

        % 3D visualisation - level 1
        function h = Plot3(obj, color)
            % plot this tube as a line as if an edge on a 3D network graph.
            % 
            % Note that datatips in plot contains specific information.
            %
            % .. todo: consider adding more info in datatips.
            %
            % Args:
            %  color: *OPTIONAL* `default = 'k'` color of line. Can be any
            %   MATLAB accepted color format, e.g. RGB.
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
            % Plot segmentation surface of tube.
            % 
            %
            % .. todo: consider adding more info in datatips.
            %
            % Args:
            %  color: *OPTIONAL* `default = 'k'` color of surface. Can be any
            %   MATLAB accepted color format, e.g. RGB.
            %  alpha (float): *OPTIONAL* `default = 0.3` opacity of surface
            %   plot.
            %  context (bool): *OPTIONAL* `default = true` also plot
            %   adjacent tubes.
            %  contextcolor: *OPTIONAL* `default = 'r'` color of adjacent 
            %   tubes. Can be any MATLAB accepted color format, e.g. RGB.
            %
            % .. todo: consider adding datatip for stats
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
                adjtubes = [obj.parent, obj.children];
                if ~isempty(obj.parent)
                    adjtubes = [adjtubes obj.parent(1).children];
                end
                for adjii = adjtubes
                    adjii.Plot3D(type=options.type, color=options.contextcolor,...
                        context=false, alpha=options.alpha)
                    hold on
                end
            end
            V = zeros(size(obj.network.skel));

            V(obj.([options.type, 'points'])) = 1;

            h = patch(isosurface(V),...
                'FaceAlpha', options.alpha,...
                'FaceColor', options.color,...
                'EdgeColor', 'none');
            obj.network.vol3daxes()


            hold off
        end

        function PlotSpline(obj,options)
            % Plot the tube spline.
            %
            % Args:
            %  color: *OPTIONAL* `default = 'k'` color of line. Can be any
            %   MATLAB accepted color format, e.g. RGB.
            %  context (bool): *OPTIONAL* `default = true` also plot
            %   adjacent tubes.
            %  contextcolor: *OPTIONAL* `default = 'r'` color of adjacent 
            %   tubes. Can be any MATLAB accepted color format, e.g. RGB.
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
            % Plot the tangental vectors of the spline at arcpoints
            % interval.
            %
            % Args:
            %  subsamp: *OPTIONAL* `default = '2'` divisive factor to
            %  subsample the arcpoints. e.g. 2 = every 2nd arcpoint tangent
            %  vector.
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
        
        function objdict = ExportSlicerLine(obj,colour)
            % Export tube object to slicer line markup format.
            %
            % .. Warning:
            %   MATLAB built in JSON encoding is sloppy with encoding
            %     floats. RBG colours as integers will not be encoded 
            %     correctly. 
            %
            
            if nargin < 2
                % default colour - turquoise
                colour = [48, 213, 200]/255;
            end


            % validate colour
            colour = validatecolor(colour);

            % get endpoints
            [X,Y,Z] = obj.I2S([obj.skelpoints(1),obj.skelpoints(end)]);
            
            % convert to global coordinates
            RAS1 = IJKtoRAS([X(1);Y(1);Z(1)]);
            RAS2 = IJKtoRAS([X(2);Y(2);Z(2)]);

            % create control point dictionaries
            cp1 = MakeControlPoint(RAS1', '1');
            cp2 = MakeControlPoint(RAS2', '2');
            
            % create measurements dict
            measurements = MakeMeasurement("length", obj.stats.euclength, "mm");
            
            % create display
            display = struct();
            display.selectedColor = colour;

            % create Line dictionary
            objdict = struct();
            objdict.type = "Line";
            objdict.coordinateSystem = "RAS";
            objdict.coordinateUnits = "mm";
            objdict.labelFormat = num2str(obj.ID);
            objdict.controlPoints = [cp1, cp2];
            objdict.locked = true;
            objdict.measurements = {measurements};
            objdict.display = display;

            % encode to json and save

            function cp = MakeControlPoint(pos, id)
                cp = struct();
                cp.id = id;
                cp.label = string(['Tube_', num2str(obj.ID),'-',id]);
                cp.associatedNodeID = "vtkMRMLScalarVolumeNode1";
                % force pos to be float
                cp.position = pos+[0.000001,0.000001,0.000001];
                % force orientation to be float
                cp.orientation = [-1.000001, 0.000001, 0.000001, 0.000001, -1.000001, 0.000001, 0.000001, 0.000001, 1.000001];
                cp.visibility = true;
            end

            function mm = MakeMeasurement(name, value, unit)
                mm = struct();
                mm.name = name;
                mm.enabled = false;
                mm.value = value;
                mm.unit = unit;
%                 mm.printFormat = " ";
            end

            function RAS = IJKtoRAS(ijk)
                % ijk must be 1 x 3 vector
                ijk = ijk + obj.network.lims(:,1);
                ijk(4) = 1;
                T = obj.network.header.Transform.T';
                RAS = T * ijk;
                RAS = RAS(1:3);
            end
        end

        function ExportOrthoPatches(obj, tarpath, casename)
            % export perpendicular slice patches of this tube.
            %
            % export the perpendicular slice patches of this tube stored in
            % source as int16 tiff files to a tar file.
            %
            %
            % .. note:
            %   This uses the GNU tar command from the system to create the
            %   final tarball. Therefore it is only expected to work on unix system.
            %   This was the only way to achieve the desired 'append'
            %   functionality that MATLAB `tar` doesn't offer.
            %
            %
            % Args:
            %   tarpath(str): path to tarfile save the exported patches.
            %       The tarfile will be be appended to or created if it 
            %       doesn't already exist.
            %   casename(char): casename to use as prefix to each patch
            %       slice name.
            % 
            % Example:
            %   >>> run CA_base.m;
            %   >>> AQnet.tubes(98).ExportOrthoPatches('patches','example')
            %
        
            % ensure char
            casename = char(casename);
            
            % choose which slices to save
            allslices = 1:length(obj.source);
            chosenslices = obj.PruneMeasure(allslices);

            % incase no slices are selected because they've all been pruned
            if isempty(chosenslices)
                warning(['No orthopatches to export. TubeID: ', num2str(obj.ID),'; Tube length: ', num2str(obj.stats.arclength),'; prune length: [', num2str(obj.prunelength'),'].'])
            end

            % save files to temp dir
            atemp_dir = tempname;
            mkdir(atemp_dir);

            % check tar path
            tarpath = parse_filename_extension(tarpath,'.tar');


            % loop through slices
            for k = chosenslices
                % save as int16
                float = obj.source{k,1};
                img = int16(float);

                % patch name
                paddedk = sprintf('%08d', k);
                patchname = [casename, ...
                    '_id_',num2str(obj.ID), ...
                    '_gen_', num2str(obj.generation), ...
                    '_slice_',paddedk];

                % save as tif to temp dir
                imgsavename = fullfile(atemp_dir, [patchname, '.tif']);
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
            
            if length(chosenslices) > 0
                % append files to tar path 
                tarcommand = ['tar --append -f ', tarpath, ' --remove-files -C ', atemp_dir,...
                    ' $(cd ', atemp_dir, ' ; echo *.tif)'];

                % execute tar
                system(tarcommand);
            end

            % remove temp dir
            rmdir(atemp_dir,'s');
        end

        function toGif(obj, filename, options)
            % Save tube patches as an animated series gif image.
            %
            % .. todo::
            %   * showcase example
            %
            % Args:
            %   filename(char): filename to save as.
            %   type(char): *OPTIONAL* `default = 'source' object property
            %       that is a 3D array to view.
            %   rings(bool array): *OPTIONAL* `default = 
            %       ones(1,size(:attr:`measures`,1))` same length as number
            %       of measures per patch, i.e. number of rings. For each
            %       ring to visualise. All rings by default.
            %   ellipses(bool): *OPTIONAL* `default = true` to view 
            %       ellipses of measure if available.
            %   points(bool): *OPTIONAL* `default = true` to view 
            %       points of measure if available.
            %   framerate(int): *OPTIONAL* `default = 20` framerate to save
            %       gif.
            %
            % Return:
            %   1 variable
            %   * s(`orthosliceViewer`): handle to created graphics
            %
            % Example:
            %   >>> run CA_base.m;
            %   >>> figure;
            %   >>> AQnet.tubes(98).toGif();
            %
            %
            %

            arguments
                obj
                filename
                options.type {mustBeMember(options.type,{'source','seg'})} = 'source'
                options.rings = ones(size(obj.measures,1),1)
                options.ellipses = true
                options.points = false
                options.framerate (1,1) = 20
            end

            % parse filename
            filename = parse_filename_extension(filename, '.gif');

            % convert from cell stack to 3D array.
            tubearray = ParseVolOut(obj, type=options.type);

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
                obj.UpdateAnnotateOrthoviewer(tubearray,ax,idx,options.rings,...
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
        
        function obj = ExportCSV(obj, filename)
            % Export tube properties and measurements to csv file.
            %
            % A list of properties from :attr:`patchprop` and measurements 
            % from :attr:`diameters` and :attr:`areas` saved into a single 
            % csv file where each slice is represented by a row. 
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
            %   >>> AQnet.tubes(98).ExportCSV('example_tube');
            %
            %

            % parse path to end in csv
            filename = parse_filename_extension(filename, '.csv');

            % properties to add - multiple measurements per slice
            tubeprops = ["diameters", "areas"];
            
            % struct properties to add - one measurement per slice
            structprops = ["patchprop"];
            
            tablestruct = struct;

            % add tube properties
            for prop = structprops
                subprops = fieldnames(obj.(prop));
                for sprop = string(subprops)'
                    tablestruct.(strcat(prop,'_',sprop)) = obj.(prop).(sprop)';
                end
            end
            
            % add tube measurements
            for prop = tubeprops
                for jj = 1:size(obj.(prop),1)
                    tablestruct.(strcat(prop,'_',string(jj))) = obj.(prop)(jj,:)';
                end
            end

            % add new row
            exporttable = struct2table(tablestruct);

            % write to csv
            writetable(exporttable, filename)

        end
        
        % utilities
        function I = S2I(obj,I1,I2,I3)
            % Fast specific implementation of MATLAB's `sub2ind
            % <https://uk.mathworks.com/help/matlab/ref/sub2ind.html>`_.
            %

            I = S2I3(size(obj.network.skel),I1,I2,I3);
        end

        function [I1,I2,I3] = I2S(obj,I)
            % Fast specific implementation of MATLAB's `ind2sub
            % <https://uk.mathworks.com/help/matlab/ref/ind2sub.html>`_.
            %

            [I1, I2, I3] = I2S3(size(obj.network.skel),I);
        end

        function obj = SetPruneLength(obj, prunelength)
            arguments
                obj
                prunelength (2,1) {mustBeNumeric}
            end
            obj.prunelength = prunelength;
        end

        function prunedprop = PruneMeasure(obj, var)
            % Prune any tube object property to the class set
            % :attr:`prunelength`.
            %
            % Uses the :attr:`patchprop`.arcpoints property to trim off
            % the desired variable property at either end.
            %
            % Args:
            %   var(matrix): n x m dimensional array where n is the number 
            %     of rings and m are the number of tube patches. 
            %
            % Return:
            %   1 variable
            %   * prunedprop(`matrix`): n x -1 dimensional array after
            %       pruning.
            %
            if isempty(obj.prunelength)
                obj.prunelength = [0 0];
            end
            nrings = size(var,1);
            pl = obj.prunelength;
            al = obj.patchprop.arcpoints;
            assert(size(al,2) == size(var,2), 'Input variable must be the same length as arclength.')

            prunebool = (al >= pl(1) & al <= (al(end) - pl(2)));
            % scale to number of measures
            prunebool = repmat(prunebool,nrings,1);
            prunedprop = var(prunebool == 1);
            prunedprop = reshape(prunedprop,nrings,[]);
        end

        function volout = ParseVolOut(obj,options)
            % Parse volume of a given property name
            %
            % See :class:`AirQuant.AirQuant`:meth:`ParseVolOut`.
            %
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
        function UpdateAnnotateOrthoviewer(obj,tubearray,ax,pos,rings,showellipses,showpoints)
            % internal function that updates the interactive plot of
            % :meth:`OrthoView`.
            %
            % delete old linetype graphics
            axesHandlesToChildObjects = findobj(ax, 'Type', 'Line');
            if ~isempty(axesHandlesToChildObjects)
                delete(axesHandlesToChildObjects);
            end
            for ii = 1:size(obj.measures,1)
                if rings(ii) == 1

                    % get centre displacement
                    canvas_sz = floor(obj.network.max_plane_sz/obj.network.plane_sample_sz);
                    image_sz = size(tubearray(:,1),1);
                    min_centre = canvas_sz/2 - image_sz/2;
                    if isa(obj.measures{ii,pos},'AQEllipse')
                        obj.measures{ii,pos}.plot(min_centre, ax, showellipses, showpoints);
                    end
                end
            end
        end
    end
end
