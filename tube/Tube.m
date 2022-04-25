% Lead Author: Ashkan Pakzad 2022. ashkanpakzad.github.io.
% See https://github.com/ashkanpakzad/AirQuant for more information.

classdef Tube < handle
    % Initialise the Tube class object.
    %
    % description
    %
    %
    % Parameters
    % ----------
    % network : :class:`AirQuant.network`
    %   tube network object that this tube is a subset.
    % parent: :class:`AirQuant.tube`
    %   parent tube objects
    % children: :class:`AirQuant.tube`
    %   child tube objects
    % localsegpoints: int, vector
    %   list of indexed points that make up the tube's local segmentation
    % skelpoints: int, vector
    %   list of indexed points that make up the tube's centreline
    % spline: struct
    %   polynomial that describes the spline
    % arclength: float, vector
    %   interval along spline in mm
    % reformedsource: float, array
    %   tube patches of source image along spline at interval specifed by arclength
    % reformedseg: float, array
    %   tube patches of seg image along spline at interval specifed by arclength
    % stats: struct
    %   stats of current tube
    %
    %
    % see also AIRWAY

    properties
        network
        relatives
        localsegpoints
        skelpoints
        spline
        arclength
        reformedsource
        reformedseg
        patchprop
        stats
    end
    properties (SetAccess = private)

    end

    methods
        function obj = Tube(network, skelpoints)
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
            % definitions
            obj.network = network;
            obj.skelpoints = skelpoints;

            obj.relatives = struct([]);
            obj.patchprop = struct([]);
            obj.stats = struct([]);
            
            obj.MakeSpline(obj)
            obj.ComputeSplinePoints(obj)
        end
        
        function obj = SetRelative(obj, tube, relation)
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
            
            obj.relatives = setfield(obj.relatives, relation, tube);
        end

        function obj = ComputeTubeCharacteristics(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add parent and child links, and calculate length of airways.
            % Save these to class to keep track of airway characteristics.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for link_index = 1:length(obj.Glink)
                % add parent idx column for each branch
                obj.Glink(link_index).parent_idx = find([obj.Glink(:).n2] == obj.Glink(link_index).n1);
                obj.Glink(link_index).child_idx = find(obj.Glink(link_index).n2 == [obj.Glink(:).n1]);
                % compute splines
                MakeSpline(obj, link_index);
                ComputeSplinePoints(obj, link_index);
            end

            % add arc and euc lengths to GLink
            [tortuosity, La, Le] = ComputeTortuosity(obj);
            La = num2cell(La);
            [obj.Glink(:).tot_arclength] = La{:};
            Le = num2cell(Le);
            [obj.Glink(:).tot_euclength] = Le{:};
            tortuosity = num2cell(tortuosity);
            [obj.Glink(:).tortuosity] = tortuosity{:};

        end

        % Classification

        function obj = SetRegion(obj, region)
            obj.RegionClassification = region;
        end

        function obj = SetGeneration(obj, gen)
            obj.GenerationClassification = gen;
        end


        % Data IO
        function obj = SaveAwyPatches(obj, prunelength)
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

        % SPLINE METHODS
        function obj = MakeSpline(obj, useparent)
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
            
            if nargin < 1
                useparent = True;
            end
            
            % get linear indexed points of previous branch if available.
            if useparent == true && ~isfield(obj.relatives,'parent')
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

        function obj = ComputeSplinePoints(obj, sample_interval)
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
            
            assert(~isempty(obj.spline), 'spline is empty, see method MakeSpline')

            if nargin < 2
                sample_interval = obj.network.spline_sample_sz;
            end
            
            % get spline points by set interval
            [obj.stats.arclength, obj.patchprop.parapoints, obj.patchprop.arcpoints] = Compute_Spline_Points(obj.spline, sample_interval);

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
            [~, CT_point_1] = AirQuant.ComputeNormal(obj.spline, ...
                obj.patchprop.parapoints(1));
            [~, CT_point_end] = AirQuant.ComputeNormal(obj.spline, ...
                obj.patchprop.parapoints(end));
            obj.stats.euclength = norm(CT_point_end - CT_point_1);
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
            obj.stats = obj.stats.arclength./obj.stats.euclength;
            assert(obj.stats >= 1, 'Impossible to get a tortuosity > 1')
        end

        % TRAVERSING AIRWAYS METHODS %%%
        % Technical: These methods traverse the airway centreline spline
        % interpolating the CT image at each spline sample point.
        function obj = CreateAirwayImage(obj, link_index)
            % Constructs perpendicular images as if travelling along an
            % airway segment in CT image and Segmentation.
            spline_t_points = ComputeSplinePoints(obj, link_index);
            % set up slice store
            TransAirwayImage = cell(length(spline_t_points),1);
            TransSegImage = cell(length(spline_t_points),1);
            for i = 1:length(spline_t_points)
                spline = obj.Splines{link_index, 1};
                % * Compute Normal Vector per spline point
                [normal, CT_point] = AirQuant.ComputeNormal(spline, ...
                    spline_t_points(i));
                % get approx airway size from distance map
                approx_diameter = ComputeDmapD(obj, CT_point);
                % compute interpolated slice size
                if obj.plane_scaling_sz ~= 0
                    plane_sz = ceil(approx_diameter*obj.plane_scaling_sz);
                else
                    plane_sz = obj.max_plane_sz;
                end
                % use max plane size if current plane size exceeds it
                if plane_sz > obj.max_plane_sz
                    plane_sz = obj.max_plane_sz;
                end
                % * Interpolate Perpendicular Slice per spline point
                [InterpAirwayImage, InterpSegImage] = ...
                    InterpolateCT(obj, normal, CT_point,  ...
                    plane_sz, obj.plane_sample_sz);
                % * Replace NaN entries in images with zero.
                InterpAirwayImage(isnan(InterpAirwayImage)) = 0;
                InterpSegImage(isnan(InterpSegImage)) = 0;
                TransAirwayImage{i,1} = InterpAirwayImage;
                TransSegImage{i,1} = InterpSegImage;
            end
            % * Save traversed image and arclength for each image
            obj.TraversedImage{link_index, 1} = TransAirwayImage;
            obj.TraversedSeg{link_index, 1} = TransSegImage;
        end


        function [CT_plane, seg_plane] = InterpolateCT(obj, normal, ...
                CT_point,  plane_sz, plane_sample_sz)
            % Interpolates a CT plane of the image.
            % Based on original function by Kin Quan 2018

            % * Construct spatial physical CT grid
            image_sz = size(obj.CT);
            [x_domain , y_domain , z_domain] = ...
                meshgrid(1:image_sz(2),1:image_sz(1),1:image_sz(3));
            x_domain = x_domain*obj.CTinfo.PixelDimensions(1);
            y_domain = y_domain*obj.CTinfo.PixelDimensions(2);
            z_domain = z_domain*obj.CTinfo.PixelDimensions(3);

            % * Get plane grid
            basis_vecs = Orthonormal_basis_with_tangent_vector(normal);
            plane_grid = Grids_coords_for_plane(basis_vecs(:,3),...
                basis_vecs(:,2), CT_point, plane_sz, plane_sample_sz);

            % * Execute cubic inperpolation on CT
            plane_intensities = interp3(x_domain,y_domain,z_domain,...
                obj.CT,plane_grid.y(:),plane_grid.x(:),...
                plane_grid.z(:),'cubic');

            % * Execute cubic inperpolation on seg
            seg_intensities = interp3(x_domain,y_domain,z_domain,...
                double(obj.seg),plane_grid.y(:),plane_grid.x(:),...
                plane_grid.z(:),'cubic');

            % Reshape
            plane_length = sqrt(length(plane_grid.y(:)));
            CT_plane = reshape(plane_intensities,...
                [plane_length plane_length]);
            seg_plane = reshape(seg_intensities,...
                [plane_length plane_length]);
        end

        function approx_diameter = ComputeDmapD(obj, CT_point)
            % Convert CT_point mm back to voxel ind
            vox_point = CT_point'./obj.CTinfo.PixelDimensions;
            % find nearest skeleton point to voxpoint
            if isempty(obj.skel_points)
                ComputeSkelPoints(obj)
            end
            P = obj.skel_points;
            k = dsearchn(P,vox_point);
            % Get radius and convert to diameter
            approx_diameter = obj.Dmap(P(k,1),P(k,2),P(k,3))*2;
            % incase of edge case, unit radius
            if approx_diameter == 0
                approx_diameter = 2;
            end
            % TODO: diameter conversion to mm
        end

        function [CT_rays, seg_rays, coords] = Raycast(obj, interpslice, interpseg, center)
            % Compute Rays
            % Getting the range limit of the ray which will be the shortest
            % distance from the centre to the bounadry Need to find the
            % limits of the raw
            image_sz = size(interpslice);
            limit_row = abs(image_sz(1) - center(1));
            limit_col = abs(image_sz(2) - center(2));
            ray_length_limit = min(limit_row, limit_col);

            %Compute rays in polar coords
            ray_angle_interval = 2*pi/obj.num_rays;
            radial = 0:obj.ray_interval:ray_length_limit;
            theata = 0:ray_angle_interval:2*pi;

            %Convert rays to cartesian coords
            x_component = radial'*cos(theata) + center(1);
            y_component = radial'*sin(theata) + center(2);

            % * Cast rays
            interpslice = double(interpslice);

            CT_rays = interp2(interpslice, x_component(:),...
                y_component(:));

            seg_rays = interp2(interpseg, x_component(:),...
                y_component(:));

            %Need to reshape
            CT_rays = ...
                reshape(CT_rays,[size(y_component,1) size(y_component,2)]);
            seg_rays = ...
                reshape(seg_rays,[size(y_component,1) size(y_component,2)]);
            coords = cat(3, x_component, y_component);
        end

        function elliptical_sturct = ComputeEllipses(obj, FWHM_points)
            %Perform the Ellipitcal fitting
            %The output will be sturct containing the major and minor lengths as well
            %as all the stop points
            elliptical_sturct = struct;
            elliptical_sturct.x_points = FWHM_points(:,1);
            elliptical_sturct.y_points = FWHM_points(:,2);
            % Compute fitting
            elliptical_sturct.elliptical_info = ...
                Elliptical_fitting(elliptical_sturct.x_points, ...
                elliptical_sturct.y_points);
            elliptical_sturct.area = ...
                elliptical_sturct.elliptical_info(3)*...
                elliptical_sturct.elliptical_info(4)*...
                pi*obj.plane_sample_sz.^2;
        end



        % LEGACY: LONG AIRWAY TAPERING METRICS.

        % long tapering
        function [logtapergrad, cum_arclength, cum_area, edgepath] = ConstructTaperPath(obj, terminal_node_idx, type)
            if nargin < 3
                type = 'inner';
            end
            switch type
                case 'inner'
                    typeidx = 1;
                case 'peak'
                    typeidx = 2;
                case 'outer'
                    typeidx = 3;
                otherwise
                    error(['Unknown type: ', type, 'please use ''inner'', ''peak'' or ''outer''.'])
            end
            % TODO: remove trachea node?
            G = digraph(obj);
            path = shortestpath(G, obj.carina_node, terminal_node_idx);

            % construct edgepath from graph2skel NB: MATLAB graph object
            % does not work.
            nodelist = [[obj.Glink(:).n1]', [obj.Glink(:).n2]'];
            edgepath = zeros(length(path)-1, 1);
            for i = 1:length(path)-1
                nodepair = [path(i), path(i+1)];
                [~,edgepath(i)] = ismember(nodepair, nodelist,'rows');
            end
            cum_arclength = 0;
            cum_area = [];
            k = 1;
            for q = 1:length(edgepath)
                i = edgepath(q);
                % skip the first of each link except the first one
                if k == 1
                    try
                        cum_area = [cum_area; obj.FWHMesl{i, typeidx}{1, 1}.area];
                        k=0;
                    catch
                        cum_area = [cum_area; NaN];
                        k=0;
                    end
                end
                % Appending arclengths
                max_current_arclength = max(cum_arclength);
                current_arclength_array = max_current_arclength+obj.arclength{i,1}(2:end);
                cum_arclength = [cum_arclength; current_arclength_array'];

                for j = 2:length(obj.arclength{i,1})
                    try
                        cum_area = [cum_area; obj.FWHMesl{i, typeidx}{j, 1}.area];
                    catch
                        cum_area = [cum_area; NaN];
                    end
                end
            end
            logtapergrad = AirQuant.ComputeTaperGrad(cum_arclength, cum_area);

            % TODO: consider using graph and edge highlight for another
            % applications.
            %h = plot(G);
            %highlight(h,'Edges',edgepath,'EdgeColor','r','LineWidth',1.5)
        end


        function terminalnodelist = ListTerminalNodes(obj)
            terminalnodelist = find([obj.Gnode(:).ep] == 1);
            % remove trachea node
            terminalnodelist(terminalnodelist ==  1) = [];
        end


        function AllTaperResults = ComputeTaperAll(obj)
            % get list of terminal branches
            terminallinklist = ListTerminalNodes(obj);
            % construct structure to save analysis results
            AllTaperResults = struct('terminalnode',num2cell(terminallinklist));
            % compute taper results for each terminal branch
            for i = 1:length(AllTaperResults)
                node_idx = AllTaperResults(i).terminalnode;

                [logtapergrad_inner, cum_arclength, cum_area_inner, edgepath] = ...
                    ConstructTaperPath(obj, node_idx, 'inner');
                [logtapergrad_peak, ~, cum_area_peak, ~] = ...
                    ConstructTaperPath(obj, node_idx, 'peak');
                [logtapergrad_outer, ~, cum_area_outer, ~] = ...
                    ConstructTaperPath(obj, node_idx, 'outer');
                AllTaperResults(i).edgepath = edgepath;
                AllTaperResults(i).arclength = cum_arclength;
                AllTaperResults(i).area_inner = cum_area_inner;
                AllTaperResults(i).area_peak = cum_area_peak;
                AllTaperResults(i).area_outer = cum_area_outer;
                AllTaperResults(i).logtapergrad_inner = logtapergrad_inner;
                AllTaperResults(i).logtapergrad_peak = logtapergrad_peak;
                AllTaperResults(i).logtapergrad_outer = logtapergrad_outer;

                try % only add lobe information if it exists
                    AllTaperResults(i).lobe = obj.Glink(obj.Gnode(node_idx).links).lobe;
                catch
                end
            end
            if isfield(AllTaperResults, 'lobe')
                AllTaperResults = sortrows(struct2table(AllTaperResults), 'lobe');
            else
                AllTaperResults = struct2table(AllTaperResults);
            end
            %             if nargout > 1
            obj.Specs.AllTaperResults = AllTaperResults;
            %             end
        end

        % SEGMENTAL ANALYSIS METHODS
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

        % TAPERING VISUALISATION METHODS
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

        % VISUALISATION
        %%% Airway Structural Tree
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

        %%% Splines
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

        % Volumetric
        function PlotSeg(obj)
            % plot segmentation
            patch(isosurface(obj.seg),'EdgeColor', 'none','FaceAlpha',0.1, 'LineStyle', 'none');
            vol3daxes(obj)
        end

        function PlotSkel(obj)
            % plot segmentation and skeleton within each other.
            patch(isosurface(obj.seg),'EdgeColor', 'none','FaceAlpha',0.1);
            hold on
            isosurface(obj.skel,'color','c')
            vol3daxes(obj)
        end

        %%% slices
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


    end
    methods (Static)
        function [normal, CT_point] = ComputeNormal(spline, point)
            % Based on original function by Kin Quan 2018
            % * interperate real point on spline
            CT_point = fnval(spline, point);
            % * get tangent of point along spline
            % differentiate along spline to get gradient
            spline_1diff = fnder(spline,1);
            tangent_vec = fnval(spline_1diff,point);
            normal = tangent_vec/norm(tangent_vec,2);
        end
    end
end
