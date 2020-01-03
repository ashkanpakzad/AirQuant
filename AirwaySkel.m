classdef AirwaySkel
    properties
        CT
        CTinfo
        seg
        branch_threshold = 0;
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
    end
    
    methods
        function obj = AirwaySkel(CTimage, CTinfo, segimage, params)
            % Initialise the AirwaySkel class.
            obj.CT = CTimage;
            obj.seg = segimage;
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
        end
        
        
        function obj = GenerateSkel(obj)
            % Generate the airway skeleton
            skel = bwskel([obj.seg]);
            [obj.Gadj,obj.Gnode,obj.Glink] =...
                Skel2Graph3D(skel,[obj.branch_threshold]);
        end
        
        
        function obj = FindTracheaCarina(obj)
            % Identify the Trachea path
            % Assumes trachea fully segmented and towards greater Z.
            [~, maxind] = max(obj.Gnode.comz);
            obj.trachea_path = obj.Gnode(maxind).links;
            obj.carina_node = obj.Gnode(maxind).conn;
        end
        
        
        function obj = CreateAirwayImage(obj, link_index)
            % Constructs perpendicular images as if travelling along an
            % airway segment.
            
            % 1. Compute whole Spline
            spline = ComputeSpline(obj, link_index);
            
            % loop along spline
            for i = 0:obj.spline_para_interval:obj.spline_para_limit;
                % 2. Compute Normal Vector per spline point
                [normal, CT_point] = ComputeNormal(spline, point);
                % 3. Interpolate Perpendicular Slice per spline point
            end
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
        
    end
    methods (Static)
        function [normal, CT_point] = ComputeNormal(spline, point)
            % Based on original function by Kin Quan 2018
            % 1. interperate real point on spline
            CT_point = fnval(spline, point);
            % 2. get tangent of point along spline
            % differentiate along spline to get gradient
            spline_1diff = fnder(spline,1);
            tangent_vec = fnval(spline_1diff,point);
            normal = tangent_vec/norm(tangent_vec,2);
        end
    end
end