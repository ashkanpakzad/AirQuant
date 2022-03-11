% Lead Author: Ashkan Pakzad 2022. ashkanpakzad.github.io.
% See https://github.com/ashkanpakzad/AirQuant for more information.

classdef TQTube < handle % handle class
    properties
        CaseDir % director for case
    end
    properties (SetAccess = private)
        SavePath % filename to load/save the object.
        %%% original Image
        CTDimRed % reduced dimensions of original CT image
        VoxDim % voxel dimensions in mm usually
        SegIndices % indices of voxels made up by segmentation
        SkelIndices % indices of voxels made up by skeletonisation
        %%% Graph Properties
        ParentNode
        ChildNode
        IgnoredLink
        LinkIndex        
        %%% Resampled image slices along graph paths/airway segments.
        TraversedImage % perpendicular CT slices of all airways
        TraversedSeg % perpendicular segmentation slices of all airways
        ArcLengths % corresponding arclength measurement traversed slices
        % Characteristics of tube
        FWHMesl % FWHMesl algorithm result for every airway {inner, peak, outer}
        Specs % Store of end metrics
        RegionClassification % e.g. lobe
        GenerationClassification % displacement from originating branch
    end
    
    methods
        %% INITIALISATION
        % Methods used to call and make AQ object.
        function obj = TQTube(LinkIndex, CaseDir, TubeType)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Initialise the TQTube class object.
            % A TQTube object, e.g. an airway or other anatomical tube.
            % Initialise object if does not exist or load if it does.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.LinkIndex = LinkIndex;
            obj.CaseDir = CaseDir;
            obj.TubeType = TubeType;
            
            name = [char(TubeType),'_',LinkIndex,'.mat'];

            obj.SavePath = fullfile(CaseDir,name);

            
            if isfile(obj.SavePath)
                disp(['Found tube stored at ', char(obj.SavePath), ' loading ...'])
                load(obj.SavePath)
                % rewrite stored path to given path
            else
                disp(['New tube to be saved at ', char(obj.SavePath), ' saving...'])
                save(obj)
            end
                
        end
        
        function obj = SetCT(obj, CTDimRed, VoxDim)
            % CTDimRed - 6 element vector of top left to bottom right
            % extremes of CT image
            % VoxDim - 3 element vector of voxel dimensions of CT in mm.
            obj.CTDimRed = CTDimRed;
            obj.VoxDim = VoxDim;
        end
        
        function obj = SetSegSkel(obj, SegIndices, SkelIndices)
            % Each is a vector of linear indices in the frame of the
            % reduced CT.
            obj.SegIndices = SegIndices;
            obj.SkelIndices = SkelIndices;
        end
        
        function obj = SetBranch(obj, ParentNode, ChildNode, IgnoredLink)
            % Set branch characteristics
            obj.ParentNode = ParentNode;
            obj.ChildNode = ChildNode;
            obj.LinkIndex = LinkIndex;
            obj.IgnoredLink = IgnoredLink;
        end
               
        function obj = ComputeBranchCharacteristics(obj)
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
                if any(link_index == obj.trachea_path)
                    continue
                end
                ComputeSpline(obj, link_index);
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
        
        %% Classification
            
        function obj = SetRegion(obj, region)
            obj.RegionClassification = region;
        end
        
        function obj = SetGeneration(obj, gen)
            obj.GenerationClassification = gen;
        end
        
            
            %% Data IO
            function save(obj)
                save(obj.SavePath, 'obj', '-v7.3')
            end
                                                
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
                
            %% HIGH LEVEL METHODS
            % methods that package lower level methods, often to apply analysis
            % to all airways rather than just individual airways.
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
            
            function obj = FindFWHMall(obj)
                % Clear all FWHMesl cells
                obj.FWHMesl = cell(length(obj.Glink),3);
                % analyse all airway segments except the trachea.
                disp('Start computing FWHM boundaries of all airway segments')
                total_branches = length(obj.Glink);
                
                % check which branches already traversed
                incomplete = cellfun(@isempty, obj.FWHMesl(:,1));
                
                for i = 1:length(obj.Glink)
                    % skip the trachea
                    if any(i == obj.trachea_path) || incomplete(i) == 0
                        continue
                    end
                    obj = FindAirwayBoundariesFWHM(obj, i);
                    disp(['FWHMesl: Completed ', num2str(i), ' of ', num2str(total_branches)])
                end
                save(obj)
                disp('FWHMesl: Done')
            end
            
            %% SPLINE METHODS
            % Technical: These methods compute the central airway spline.
            function ComputeSpline(obj, link_index)
                % Computes a smooth spline of a single graph edge.
                % Based on original function by Kin Quan 2018
                
                % The input is the list of ordered index
                % The output is the smooth spline as a matlab sturct
                
                % get linear indexed points of current and previous branch,
                % combine.
                previous_awy = find([obj.Glink(:).n2] == obj.Glink(link_index).n1);
                [x_p1, y_p1, z_p1] = ind2sub(size(obj.CT), obj.Glink(previous_awy).point);
                [x_p2, y_p2, z_p2] = ind2sub(size(obj.CT), obj.Glink(link_index).point);
                x_point = [x_p1, x_p2];
                y_point = [y_p1, y_p2];
                z_point = [z_p1, z_p2];
                
                %Smooth all points using moving average
                voxel_sz = obj.CTinfo.PixelDimensions;
                smooth_x = smooth(x_point*voxel_sz(1),11, 'moving');
                smooth_y = smooth(y_point*voxel_sz(2),11, 'moving');
                smooth_z = smooth(z_point*voxel_sz(3),11, 'moving');
                
                % extract just current airway smoothed points
                csmooth_x = smooth_x(length(x_p1)+1:end);
                csmooth_y = smooth_y(length(x_p1)+1:end);
                csmooth_z = smooth_z(length(x_p1)+1:end);
                
                %Complete smooth data
                smooth_data_points = [csmooth_x csmooth_y csmooth_z]';
                
                %Generating the spline
                obj.Splines{link_index, 1} = cscvn(smooth_data_points);
            end
            
            function t_points = ComputeSplinePoints(obj, link_index)
                % * Compute Spline if necessary
                if isempty(obj.Splines{link_index, 1})
                    ComputeSpline(obj, link_index)
                end
                spline = obj.Splines{link_index, 1};
                
                [t_points, arc_length] = spline_points(spline, obj.spline_sample_sz);
                
                % save parametrised sample points in second column.
                obj.Splines{link_index, 2} = t_points;
                
                % save arc_length measurement of each sample point.
                obj.arclength{link_index, 1} = arc_length;
                
            end
            
            %% TRAVERSING AIRWAYS METHODS %%%
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
            
            %% MEASUREMENT METHODS
            % Methods that measure the airway size on interpolated CT images.
            function obj = FindAirwayBoundariesFWHM(obj, link_index, outlier_removal)
                if nargin < 3
                    outlier_removal = true;
                end
                %Based on function by Kin Quan 2018 that is based on Kiraly06
                
                slices_sz = size(obj.TraversedImage{link_index, 1}, 1);
                
                % Prepping the outputs
                raycast_FWHMl = cell(slices_sz, 1);
                raycast_FWHMp = cell(slices_sz, 1);
                raycast_FWHMr = cell(slices_sz, 1);
                raycast_center = cell(slices_sz, 1);
                
                
                % For every traversed slice
                for k = 1:slices_sz
                    try
                        % * Compute airway centre
                        center =  AirwayCenter(...
                            obj.TraversedImage{link_index, 1}{k,1}, ...
                            obj.TraversedSeg{link_index, 1}{k,1});
                        
                        % * Raycast
                        [CT_rays, seg_rays, coords]= Raycast(obj, ...
                            obj.TraversedImage{link_index, 1}{k,1}, ...
                            obj.TraversedSeg{link_index, 1}{k,1}, center);
                        
                        % * Compute FWHM
                        [FWHMl, FWHMp, FWHMr] = AirQuant.computeFWHM(CT_rays,...
                            seg_rays, coords, outlier_removal);
                        
                        % * Compute Ellipses
                        FWHMl_ellipse = ComputeEllipses(obj, FWHMl);
                        FWHMp_ellipse = ComputeEllipses(obj, FWHMp);
                        FWHMr_ellipse = ComputeEllipses(obj, FWHMr);
                        
                        % * Record, catch incase a slice fails.
                        raycast_FWHMl{k,1} = FWHMl_ellipse;
                        raycast_FWHMp{k,1} = FWHMp_ellipse;
                        raycast_FWHMr{k,1} = FWHMr_ellipse;
                        raycast_center{k,1} = center;
                    catch
                        % segmentation exceeds interpolated slice therefore no
                        % measurement recorded.
                        raycast_FWHMl{k,1} = NaN;
                        raycast_FWHMp{k,1} = NaN;
                        raycast_FWHMr{k,1} = NaN;
                        raycast_center{k,1} = NaN;
                    end
                end
                obj.FWHMesl{link_index, 1} = raycast_FWHMl;
                obj.FWHMesl{link_index, 2} = raycast_FWHMp;
                obj.FWHMesl{link_index, 3} = raycast_FWHMr;
                obj.FWHMesl{link_index, 4} = raycast_center;
                
            end
            
            function [CT_rays, seg_rays, coords] = Raycast(obj, interpslice, interpseg, center)
                % * Compute Rays
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
            
            %% EXTENT Methods
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
            
            function [overextentidx, compare] = compareextent(obj, populationstats, threshold)
                % threshold = minimum gen expected
                % populationstats = table generated by summarising results from
                % genperlobe over several images.
                subject_genperlobe = table2array(extentstats(obj));
                if istable(populationstats)
                    populationstats = table2array(populationstats);
                end
                gendiff = size(populationstats, 1) - size(subject_genperlobe, 1);
                if gendiff < 0
                    populationstats = [populationstats; zeros(abs(gendiff),6)];
                elseif gendiff > 0
                    subject_genperlobe = [subject_genperlobe; zeros(abs(gendiff),6)];
                end
                
                % compare this subject to population
                compare = floor(subject_genperlobe - populationstats);
                
                overextentidx = false(length(obj.Glink),1);
                %             safenodes = [];
                varNames = AirQuant.lobelabels();
                % BF search
                [~,bfs_Gind] = bfsearch(obj.Gdigraph, 1,'edgetonew');
                % convert graph idx into link index
                bfs_Lind = zeros(size(bfs_Gind));
                for ind = 1:length(bfs_Gind)
                    bfs_Lind(ind) = obj.Gdigraph.Edges.Label(bfs_Gind(ind));
                end
                % reverse order
                bfs_Lind = flip(bfs_Lind);
                for lobe = 1:6
                    lobeid = varNames{1,lobe};
                    for gen = size(subject_genperlobe, 1):-1:threshold+1
                        n_overextent = compare(gen, lobe);
                        if n_overextent <= 0
                            continue
                        end
                        overextentawy(obj, gen, lobeid, n_overextent);
                    end
                end
                %             for awy = 1:length(obj.Glink)
                %                 if ~isempty(ismember(obj.Glink(awy).generation,[1:threshold]))
                %                     overextentidx(awy) = false;
                %                 end
                %             end
                %             function overextentawy(obj, gen, lobe, n_overextent)
                %             Nawy = length(obj.Glink);
                %             counts = 1;
                %             mode = 0;
                %             for i = 1:Nawy
                %                 if counts > n_overextent
                %                     mode = 1; % safe mode
                %                 end
                %                 checksafe = isempty(ismember([obj.Glink(i).n1, obj.Glink(i).n2], safenodes));
                %                 if obj.Glink(i).generation == gen && strcmp(obj.Glink(i).lobe, lobe) && checksafe == 0
                %                     if mode == 0
                %                         overextentidx(i) = true; % highlight as excess
                %                         counts = counts + 1;
                %                     elseif mode == 1
                %                         safenodes = [safenodes; obj.Glink(i).n1; obj.Glink(i).n2];
                %                     end
                %                 end
                %             end
                %             end
                %         end
                
                function overextentawy(obj, gen, lobe, n_overextent)
                    counts = 1;
                    for i = 1:length(bfs_Lind)
                        if counts > n_overextent
                            break
                        end
                        awyind = bfs_Lind(i);
                        if obj.Glink(awyind).generation == gen && strcmp(obj.Glink(awyind).lobe, lobe)
                            overextentidx(awyind) = true; % highlight as excess
                            counts = counts + 1;
                        end
                    end
                end
            end
            
            function h = plotoverextent(obj, overextentidx, type)
                overextentedge = overextentidx(obj.Gdigraph.Edges.Label);
                h = plot(obj, type);
                highlight(h,'Edges',overextentedge, 'LineStyle', ':', 'LineWidth',1)
            end
            
            %% LEGACY: LONG AIRWAY TAPERING METRICS.
            
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
            
            %% SEGMENTAL ANALYSIS METHODS
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
            
            function lobar_intertaper = ComputeLobarInterTaper(obj, prunelength)
                % TODO: great potential for improving efficiency.
                % Compare every lobar airway with the diameter of the first
                % airway of that lobe.
                
                if nargin < 2
                    prunelength = [0 0];
                end
                
                % Check if lobe algorithm was successful, fail if not.
                if ~isfield(obj.Glink, 'lobe')
                    warning('ComputeLobarInterTaper failed, run lobe classification successfully first.')
                end
                
                % use output from intrataperall
                [~, averagediameter] = ComputeIntraTaperAll(obj, prunelength);
                indexed_lobes = convertCharsToStrings([obj.Glink.lobe]);
                % loop through branches
                labels = AirQuant.LobeLabels();
                lobar_intertaper = NaN(length(obj.arclength), 3);
                for iii = 1:length(labels)
                    current_lobe = labels(iii);
                    lobe_idx = find(indexed_lobes == current_lobe);
                    
                    for ii = lobe_idx
                        if obj.Glink(ii).generation < 3
                            continue
                        end
                        
                        % identify avg vol of 2nd lobe gen
                        s2nd_branch = find([obj.Glink.generation] == 2 & indexed_lobes == current_lobe);
                        for jj = 1:3
                            s2nd_vol = mean(averagediameter(s2nd_branch, jj),'omitnan');
                            % identify parent by predecessor node
                            lobar_intertaper(ii,jj) = ((s2nd_vol - averagediameter(ii,jj))...
                                /(s2nd_vol)) * 100;
                        end
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
            
            %% TAPERING VISUALISATION METHODS
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
            
            %% Graph Methods
            function G = digraph(obj)
                % TODO: REDUNDANT FUNCTION?
                % compute edge weights
                weights = zeros(length(obj.Glink), 1);
                for i = 1:length(obj.Glink)
                    weights(i) = length(obj.Glink(i).point);
                end
                % add edges
                G = digraph([obj.Glink(:).n1],[obj.Glink(:).n2], weights);
            end
            
            %% VISUALISATION
            %%% Airway Strucutral Tree
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

            %%% Splines
            function PlotSplineTree(obj)
                % loop through every branch, check spline has already been
                % computed, compute if necessary. Skip trachea. Plot spline.
                % TODO: This may not work as well when VoxelSize ~= [1,1,1]
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
                        % TODO: rewrite this bit....
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
            %% EXPORT METHODS
            % methods for exporting processed data.
            function exportCT(obj, savename)
                if nargin == 1
                    savename = 'AQCT_export';
                end

                % reduce datatype and change header info
                AQ_CT = int16(obj.CT);
                header = obj.CTinfo;
                header.Description = 'CT exported from AirQuant';
                header.ImageSize = size(AQ_CT);
                niftiwrite(AQ_CT, savename, header, 'Compressed', true);
            end
            
            function exportlobes(obj, savename)
                if nargin == 1
                    savename = 'AQlobes_export';
                end
                % export airway segmentation labelled by lobes to nii.gz
                
                % classify by branch first
                branch_seg = ClassifySegmentation(obj);
                lobeid = {'B','RUL','RML','RLL','LUL','LML','LLL'};
                lobe_airway_seg = zeros(size(obj.seg));
                
                % convert branch classification to lobe classification
                try
                    for i = 1:length(obj.Glink)
                        lobe_airway_seg(branch_seg == i) = find(strcmp(lobeid, obj.Glink(i).lobe));
                    end
                catch
                    error('Need to run ComputeAirwayLobes first')
                end
                
                % reduce datatype and change header info
                lobe_airway_seg = uint8(lobe_airway_seg);
                header = obj.CTinfo;
                header.Datatype = 'uint8';
                header.BitsPerPixel = '8';
                header.Description = 'Airway Lobe Segmentation using AirQuant, layers: B,RUL,RML,RLL,LUL,LML,LLL';
                header.ImageSize = size(lobe_airway_seg);

                niftiwrite(lobe_airway_seg, savename, header, 'Compressed', true);
            end
            
            function innerareas = GetAreas(obj)
                % Extract the inner area of all branches from FWHM results and output
                
                if  ~isfield(obj.Specs, 'innerareas')
                    innerareas = cell(length(obj.Glink), 1);
                    for ii = 1:length(obj.Glink)
                        for jj = 1:length(obj.FWHMesl{ii,1})
                            try
                                innerareas{ii, 1} = [innerareas{ii, 1}, obj.FWHMesl{ii,1}{jj, 1}.area];
                            catch
                                innerareas{ii, 1} = [innerareas{ii, 1}, NaN];
                            end
                        end
                    end
                    maxlen = max(cellfun(@length, innerareas));
                    innerareas = cellfun(@(x)([x nan(1, maxlen - length(x))]), innerareas, 'UniformOutput', false);
                    innerareas = cell2mat(innerareas);
                    obj.Specs.innerareas = innerareas;
                else
                    innerareas = obj.Specs.innerareas;
                end
                
            end
            
            function innerdiameters = GetDiameters(obj)
                % Extract the inner area of all branches from FWHM results,
                % convert to diameter and output.
                if  ~isfield(obj.Specs, 'innerdiameters')
                    innerareas = GetAreas(obj);
                    innerdiameters = sqrt(innerareas/pi)*2;
                    obj.Specs.innerdiameters = innerdiameters;
                else
                    innerdiameters = obj.Specs.innerdiameters;
                end
            end
            
            function exportraw(obj, savename)
                % export Glink table and arclengths and inner diameter
                % measurements
                
                % Glink to table
                BranchInfo = struct2table(obj.Glink);
                
                % arclengths to table
                arclengths = obj.arclength;
                maxlen = max(cellfun(@length, arclengths));
                arclengths = cellfun(@(x)([x nan(1, maxlen - length(x))]), arclengths, 'UniformOutput', false);
                arclengths = cell2mat(arclengths);
                
                % diameters
                innerdiameters = GetDiameters(obj);
                
                % export csvs
                writetable(BranchInfo, ['BranchInfo_',savename, '.csv']);
                writematrix(arclengths, ['ArcLengths_',savename, '.csv']);
                writematrix(innerdiameters, ['InnerDiameters_',savename, '.csv']);
                
            end
        end
        %% STATIC METHODS
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
            
            
            function [logtapergradient, displacement] = ComputeTaperGrad(arclength, area)
                % identify NaN data points
                idx = isnan(area);
                % compute logtapergrad, ignoring nan values
                p_coeff = polyfit(arclength(~idx),log(area(~idx)),1);
                % positive if thinning, therefore multiply by -1.
                logtapergradient = p_coeff(1) * -1; % logtapergradient
                if nargout > 1
                    displacement = p_coeff(2); % displacement
                end
            end
            
            
            function [FWHMl, FWHMp, FWHMr] = computeFWHM(CT_rays, seg_rays, coords, outlier_removal)
                if nargin < 4
                    outlier_removal = true;
                end
                %This is to perfrom the ray casting measurents - the input is in a sturct
                % as
                
                % The basis of the code is base on the description - Virtual Bronchoscopy for Quantitative Airway Analysis
                % by A P Kiraly et al.
                
                %Need to find the length of the ray - this need a loop
                number_of_rays = size(CT_rays,2);
                
                FWHMl_r = nan(number_of_rays, 1);
                FWHMp_r = nan(number_of_rays, 1);
                FWHMr_r = nan(number_of_rays, 1);
                
                %performing a loop
                for ray = 1:number_of_rays
                    
                    CT_profile = CT_rays(:,ray);
                    seg_profile = seg_rays(:,ray);
                    
                    % * Find the edge of the interpolated segmentation.
                    if seg_profile(1) < 0.5
                        continue % skip if seg edge does not exist
                    end
                    
                    % Identify the last vaule that is above the 0.5
                    ind_ray = (seg_profile < 0.5);
                    seg_half = find(ind_ray,1);
                    
                    if all(ind_ray ~= 1)
                        continue % skip if fail to find point <0.5
                    end
                    
                    % * Find FWHM peak
                    [max_int_array , max_location_array] = ...
                        findpeaks(CT_profile);
                    
                    if isempty(max_int_array)
                        continue % skip ray if peak does not exist
                    end
                    
                    % identify the peak closest to seg half point
                    index_diff = abs(max_location_array - seg_half);
                    [~,closest_max] = min(index_diff);
                    % ensure the point is unique
                    FWHMp = max_location_array(closest_max(1));
                    
                    % * Compute FWHM on left side
                    % loop through profile from FWHM peak to centre
                    % find first minima from right to left.
                    for i = FWHMp:-1:2
                        if ~(CT_profile(i) >= CT_profile(i - 1))
                            break
                        end
                        i = 1; % incase inner is at beginning of profile
                    end
                    FWHMi = i; % inner FWHM curve
                    
                    threshold_int_left = ...
                        (CT_profile(FWHMp) + CT_profile(FWHMi))/2;
                    FWHMl = Finding_midpoint_stop(CT_profile,threshold_int_left,FWHMi,FWHMp);
                    
                    % * Compute FWHM on right side
                    % loop through profile from FWHM peak to distal
                    % find first minima from left to right.
                    for i = FWHMp:length(CT_profile)-1
                        if CT_profile(i) <= CT_profile(i + 1)
                            break
                        end
                        i = length(CT_profile);
                    end
                    FWHMo = i;  % outer FWHM curve
                    % TODO: move threshold calculation into function.
                    threshold_int_right = ...
                        (CT_profile(FWHMp) + CT_profile(FWHMo))/2;
                    FWHMr = Finding_midpoint_stop_right(CT_profile,threshold_int_right,FWHMo,FWHMp);
                    
                    % concat points together
                    FWHMl_r(ray) = FWHMl;
                    FWHMp_r(ray) = FWHMp;
                    if isempty(FWHMr) % in case right peak not within profile
                        FWHMr = NaN;
                    end
                    FWHMr_r(ray) = FWHMr;
                    
                end
                
                % anomaly correction on inner ONLY
                valid_rays = [1:number_of_rays]';
                if outlier_removal == true
                    [~,anomalies] = rmoutliers(FWHMl_r, 'median');
                    valid_rays = valid_rays(~anomalies);
                    FWHMl_r = FWHMl_r(~anomalies);
                    FWHMp_r = FWHMp_r(~anomalies);
                    FWHMr_r = FWHMr_r(~anomalies);
                end
                
                % convert radial index into plane coordinates
                FWHMl_x = DualIndex(coords(:,:,1), FWHMl_r, valid_rays);
                FWHMl_y = DualIndex(coords(:,:,2), FWHMl_r, valid_rays);
                FWHMp_x = DualIndex(coords(:,:,1), FWHMp_r, valid_rays);
                FWHMp_y = DualIndex(coords(:,:,2), FWHMp_r, valid_rays);
                FWHMr_x = DualIndex(coords(:,:,1), FWHMr_r, valid_rays);
                FWHMr_y = DualIndex(coords(:,:,2), FWHMr_r, valid_rays);
                
                % output
                FWHMl = cat(2, FWHMl_x, FWHMl_y);
                FWHMp = cat(2, FWHMp_x, FWHMp_y);
                FWHMr = cat(2, FWHMr_x, FWHMr_y);
            end
            
            function labels = LobeLabels()
                % returns cell of lobe lables used in AirQuant
                labels = {'RUL','RML','RLL','LUL','LML','LLL'};
            end
        end
end