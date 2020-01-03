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
          % Compute Spline
          % Compute Normal Vectors
          % Interpolate Perpendicular Slices
      end
      
      
   end
end