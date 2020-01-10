classdef AirwaySkel
    properties
        CT
        CTinfo
        seg
        branch_threshold = 3;
        % CT Properties/resampling params
        physical_plane_length = 40;
        physical_sampling_interval = 0.3;
        spline_sampling_interval = 0.25;
        % Ray params
        rays = 50;
        ray_interval = 0.2;
        nan_replace_int = 0;
    end
    properties (SetAccess = private)
        % Graph Properties
        Gadj
        Gnode
        Glink
        trachea_path
        carina_node
        % Resampled image slices along graph paths/airway segments.
        TraversedImage
        TraversedSeg
        arclength
        
    end
    
    methods
        function obj = AirwaySkel(CTimage, CTinfo, segimage, params)
            % Initialise the AirwaySkel class.
            obj.CT = CTimage;
            % TODO: consider preprocess segmentation to keep largest
            % connected component.
            % ensure no holes in segmentation
            obj.seg = imfill(segimage,'holes');
            obj.CTinfo = CTinfo;
            % set params
            if ~isempty(params)
                obj.branch_threshold = params.branch_threshold;
                obj.physical_plane_length = params.physical_plane_length;
                obj.physical_sampling_interval = params.physical_sampling_interval;
                obj.spline_sampling_interval = params.spline_sampling_interval;
                obj.rays = params.rays;
                obj.ray_interval = params.ray_interval;
                obj.nan_replace_int = params.nan_replace_int;
            end
            % graph airway skeleton
            obj = GenerateSkel(obj);
            % identify trachea and carina
            obj = FindTracheaCarina(obj);
            % set up empty cell for traversed images
            obj.TraversedImage = cell(length(obj.Glink),1);
            % set up empty specs doubles
            obj.arclength = cell(length(obj.Glink),1);
            %obj.specs(length(obj.Glink)) = struct();
        end
        
        
        function obj = GenerateSkel(obj)
            % Generate the airway skeleton
            skel = bwskel(obj.seg,'MinBranchLength', obj.branch_threshold);
            [obj.Gadj,obj.Gnode,obj.Glink] =...
                Skel2Graph3D(skel,0);
        end
        
        
        function obj = FindTracheaCarina(obj)
            % Identify the Trachea path and the carina node.
            % Assumes trachea fully segmented and towards greater Z.
            
            % Smoothen
            [~, maxind] = max([obj.Gnode.comz]);
            obj.trachea_path = obj.Gnode(maxind).links;
            obj.carina_node = obj.Gnode(maxind).conn;
        end
        
        
        function obj = AirwayImageAll(obj)
            % Traverse all airway segments except the trachea.
            disp('Start traversing all airway segments')
            total_branches = length(obj.Glink);
            for i = 1:length(obj.Glink)
                % skip the trachea
                if i == obj.trachea_path
                    continue
                end
                obj = CreateAirwayImage(obj, i);
                disp(['Completed ', num2str(i), ' of ', num2str(total_branches)])
            end
        end
        
        %%%
        function obj = CreateAirwayImage(obj, link_index)
            % Constructs perpendicular images as if travelling along an
            % airway segment in CT image and Segmentation.
            
            % * Compute whole Spline
            spline = ComputeSpline(obj, link_index);
            spline_para_limit = spline.breaks(end);
            spline_points = 0:obj.spline_sampling_interval:spline_para_limit;
            
            % loop along spline
            % TODO: auto calc 133
            TransAirwayImage = zeros(133,133,length(spline_points));
            TransSegImage = zeros(133,133,length(spline_points));
            arc_length = zeros(length(spline_points),1);
            for i = 1:length(spline_points)
                % * Compute Normal Vector per spline point
                [normal, CT_point] = AirwaySkel.ComputeNormal(spline, spline_points(i));
                % * Interpolate Perpendicular Slice per spline point
                [TransAirwayImage(:,:,i), TransSegImage(:,:,i)] = ...
                    InterpolateCT(obj, normal, CT_point);
                % * Compute real arc_length at this spline point
                arc_length(i) = Arc_length_to_point(spline_points(i),spline);
            end
            % * Replace NaN entries in images with zero.
            TransAirwayImage(isnan(TransAirwayImage)) = 0;
            TransSegImage(isnan(TransAirwayImage)) = 0;
            % * Save traversed image and arclength for each image
            obj.TraversedImage{link_index, 1} = TransAirwayImage;
            obj.TraversedSeg{link_index, 1} = TransSegImage;
            obj.arclength{link_index, 1} = arc_length;
        end
        
        
        function [CT_plane, seg_plane] = InterpolateCT(obj, normal, CT_point)
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
                basis_vecs(:,2), CT_point, obj.physical_plane_length,...
                obj.physical_sampling_interval);
            
            % * Execute cubic inperpolation on CT
            plane_intensities = interp3(x_domain,y_domain,z_domain,...
                obj.CT,plane_grid.y(:),plane_grid.x(:),...
                plane_grid.z(:),'cubic');
            
            % * Execute cubic inperpolation on CT
            seg_intensities = interp3(x_domain,y_domain,z_domain,...
                double(obj.seg),plane_grid.y(:),plane_grid.x(:),...
                plane_grid.z(:),'cubic');
            
            % Reshape
            % TODO: Look at what these two lines do.
            plane_length = sqrt(length(plane_grid.y(:)));
            CT_plane = reshape(plane_intensities,...
                [plane_length plane_length]);
            seg_plane = reshape(seg_intensities,...
                [plane_length plane_length]);
        end
        
        
        function spline = ComputeSpline(obj, link_index)
            % Computes a smooth spline of a single graph edge.
            % Based on original function by Kin Quan 2018
            
            %The input is the list of ordered index
            %The output is the smooth spline as a matlab sturct
            
            %Convert into 3d Points
            [x_point, y_point, z_point] = ind2sub(size(obj.CT),obj.Glink(link_index).point);
            
            %Smoothing the data
            voxel_size = obj.CTinfo.PixelDimensions;
            smooth_x = smooth(x_point*voxel_size(1));
            smooth_y = smooth(y_point*voxel_size(2));
            smooth_z = smooth(z_point*voxel_size(3));
            
            %Complete smooth data
            smooth_data_points = [smooth_x smooth_y smooth_z]';
            
            %Generating the spline
            spline = cscvn(smooth_data_points);
        end
        
        %%% TAPERING MEASUREMENTS %%%
        FindAirwayBoundariesFWHM(obj, link_index)
        %Based on function by Kin Quan 2018 that is based on Kiraly06 
        
        number_of_slice = size(tapering_image,3);
        
        %Getting the pixel size - this is to convert it into the phyical diamter
        pixel_size = plane_input_sturct.physical_sampling_interval.^2;
        
        %Getting the outputs
        ellptical_info_cell = {};
        ellptical_info_cell_wall = {};
        ellptical_area = [];
        ellptical_area_wall = [];
        slice_fail_case = [];
        slice_fail_case_wall = [];
        
        %% Perfroming the loop
        for k = 1:size(obj.TraversedImage{link_index, 1}, 3)
            
            % * Check that airway centre is slice centre
            
            % Recompute centre if necessary
            
            
            
            
            plane_input_sturct.cross_sectional_image = binary_tapering_seg(:,:,k);
            plane_input_sturct.center = ...
                Return_centre_pt_image(plane_input_sturct.cross_sectional_image);
            %PLacing the raw image inorder to give a more
            plane_input_sturct.raw_sectional_image = tapering_image(:,:,k);
            
            %% This is to look at the centre of the segmentation
            [centre_ind , new_centre] =  ...
                Check_centre_with_segmentation(plane_input_sturct.cross_sectional_image,...
                binary_tapering_seg(:,:,k));
            
            if ~centre_ind
                plane_input_sturct.center = fliplr(new_centre);
            end
            
            % Raycast to find the lumen boundary
            
            %Computing the area with the describeed method
            [area_info_sturct_wall, area_info_sturct] =...
                Cross_sectional_FWHM_SL(plane_input_sturct);
            try
                %Recording the infomation
                ellptical_info_cell = ...
                    cat(1,ellptical_info_cell,area_info_sturct);
                
                %Recording the area in phyical information
                single_area = (area_info_sturct.area)*pixel_size;
                ellptical_area = cat(1,ellptical_area,single_area);
                
            catch
                %This is the case when the ray casting fails
                area_info_sturct = NaN;
                ellptical_info_cell = ...
                    cat(1,ellptical_info_cell,area_info_sturct);
                
                %Recording the area in phyical information
                single_area = NaN;
                ellptical_area = cat(1,ellptical_area,single_area);
                
                %Recroding the failed case
                slice_fail_case = cat(1,slice_fail_case,k);
                
            end
            %%  safety catch - wall
            try
                
                %Recording the infomation
                ellptical_info_cell_wall = ...
                    cat(1,ellptical_info_cell_wall,area_info_sturct_wall);
                
                single_area_wall = (area_info_sturct_wall.area)*pixel_size;
                ellptical_area_wall = cat(1,ellptical_area_wall,single_area_wall);
                
            catch
                %This is the case when the ray casting fails
                area_info_sturct_wall = NaN;
                ellptical_info_cell_wall = ...
                    cat(1,ellptical_info_cell_wall,area_info_sturct_wall);
                
                single_area_wall = NaN;
                ellptical_area_wall = cat(1,ellptical_area_wall,single_area_wall);
                
                %Recroding the failed case
                slice_fail_case_wall = cat(1,slice_fail_case_wall,k);
            end
            
        end
        
        %% Getting the outputs
        tapering_info_sturct = struct;
        tapering_info_sturct.elliptical_info = ellptical_info_cell;
        tapering_info_sturct.elliptical_info_wall = ellptical_info_cell_wall;
        tapering_info_sturct.phyiscal_area = ellptical_area;
        tapering_info_sturct.phyiscal_area_wall = ellptical_area_wall;
        tapering_info_sturct.fail_cases = slice_fail_case;
        tapering_info_sturct.fail_cases_wall = slice_fail_case_wall;
        
        
        
        %%% VISUALISATION %%%
        function PlotTree(obj)
            % Plot the airway tree with nodes and links
            % Original Function by Ashkan Pakzad on 27th July 2019.
            
            X = [obj.Gnode.comy];
            Y = [obj.Gnode.comx];
            Z = [obj.Gnode.comz];
            nums = string(1:length(X));
            
            isosurface(bwskel([obj.seg]));
            hold on
            plot3(X,Y,Z, 'r.', 'MarkerSize', 15);
            text(X+1,Y+1,Z+1, nums)
            axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            view(80,0)
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