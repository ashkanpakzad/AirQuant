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
        ID
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
        stats
        savename
    end
    properties (SetAccess = private)

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

        function obj = SetParent(obj,tube)
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

        % classification
        function obj = SetRegion(obj, region)
            obj.RegionClassification = region;
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

        % visualisation - tapering
        function PlotMeasure(obj, terminal_node_idx, type)
            % Set relative to current tube object.
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

        %%% visualisation - slices
        function PlotAirway3(obj, link_index)
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
                %                         obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(2)+min_centre,'y');
                %TODO: set third colour more appropiately

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

        function s = OrthoView(obj, options)
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
            end

            % convert from cell stack to 3D array.
            tubearray = tubestack(obj, type=options.type);
            % display with orthoview
            fig = figure;
            s = orthosliceViewer(tubearray, 'DisplayRangeInteraction','off', ...
                'ScaleFactors',[obj.network.plane_sample_sz, ...
                obj.network.plane_sample_sz, ...
                obj.network.spline_sample_sz], 'CrosshairLineWidth', 0.3);
            % Can only alter size of figure window after orthosliceviewer.
            fig.Name = 'AirQuant: Ortho View';
            fig.Units = 'normalized';
            fig.Position = [0.1,0.01,0.6,0.9];
        end

        function h = ReformatAirway(obj,slice_idx)
            % short desc
            %
            % long desc
            %
            % .. todo: add documentation to this function
            %       * Needs attention
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            % get reformatted airway stack
            tubearray = tubestack(obj);
            if nargin < 3
                % set default to middle
                slice_idx = round(size(tubearray,1)/2);
            end
            % generate image
            img = squeeze(tubearray(slice_idx,:,:));
            x = [0 obj.arclength{link_index,1}(end)];
            y = [0 obj.max_plane_sz];
            h = imagesc(x, y, img);
            colormap('gray')
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
        
        function obj = save(obj)
            % assign unique serialnumber to object if not existing
            % delete network
            % save object within network folder
            % readd network 
        end

        function obj = load(obj)
            % load object
            % attach to network
            % add tubes to object
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

        function tubearray = tubestack(obj,options)
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
end
