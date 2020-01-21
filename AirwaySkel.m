classdef AirwaySkel
    properties
        CT
        CTinfo
        seg
        branch_threshold = 2;
        % CT Properties/resampling params
        physical_plane_length = 40;
        physical_sampling_interval = 0.3;
        spline_sampling_interval = 0.25;
        % Ray params
        num_rays = 50;
        ray_interval = 0.2;
    end
    properties (SetAccess = private)
        % Graph Properties
        Gadj
        Gnode
        Glink
        Gdigraph
        trachea_path
        carina_node
        % Resampled image slices along graph paths/airway segments.
        TraversedImage
        TraversedSeg
        arclength
        FWHMesl
        Specs
    end
    
    methods
%% INITIALISATION METHODS
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
                obj.num_rays = params.num_rays;
                obj.ray_interval = params.ray_interval;
            end
            % graph airway skeleton
            obj = GenerateSkel(obj);
            % identify trachea and carina
            obj = FindTracheaCarina(obj);
            % Convert into digraph
            obj = AirwayDigraph(obj);
            % set up empty cell for traversed images
            obj.TraversedImage = cell(length(obj.Glink),1);
            % set up empty specs doubles
            obj.arclength = cell(length(obj.Glink),1);
            obj.Specs = struct();
            % set up empty cell for recording raycast/fwhmesl method
            obj.FWHMesl = cell(length(obj.Glink),3);
        end
        
        
        function obj = GenerateSkel(obj)
            % Generate the airway skeleton
            skel = bwskel(obj.seg,'MinBranchLength', obj.branch_threshold);
            [obj.Gadj,obj.Gnode,obj.Glink] =...
                Skel2Graph3D(skel,0);

        end
        
        
        function obj = FindTracheaCarina(obj)
            % Identify the Trachea path and the carina node.FindFWHMall
            % Assumes trachea fully segmented and towards greater Z.
            
            % Smoothen
            [~, maxind] = max([obj.Gnode.comz]);
            obj.trachea_path = obj.Gnode(maxind).links;
            obj.carina_node = obj.Gnode(maxind).conn;
        end
        
        
        function obj = AirwayDigraph(obj)
            % To convert the Airway graph into a digraph from carina to
            % distal.
            
            % TODO: consider identifying trachea node and working out
            % carina by connectivity of the digraph.
            
            % Create digraph with edges in both directions, loop through
            % and remove if found later in the BF search.
            G = digraph(obj.Gadj);
            % BF search from carina node to distal.
            node_discovery = bfsearch(G,obj.carina_node);
            % half of the edges will be removed
            removal = zeros(height(G.Edges)/2 , 1);
            j = 1;
            for i = 1:height(G.Edges)
                % identify the rank in order of discovery
                rank1 = find(~(node_discovery-G.Edges.EndNodes(i,1)));
                rank2 = find(~(node_discovery-G.Edges.EndNodes(i,2)));
                if rank1 > rank2
                    % record if not central-->distal
                    removal(j) = i;
                    j = j + 1;
                end      
            end
            G = rmedge(G,removal);
            obj.Gdigraph = G;
            
            % Need to ensure directions in Glink are also in the correct
            % direction. Swap over if not.
            for i = 1:length(obj.Glink)
                glink_nodes = [obj.Glink(i).n1, obj.Glink(i).n2];
                % if not central-->distal
                if ~ismember(glink_nodes, G.Edges.EndNodes, 'rows')
                    % Swap round node assignments
                    obj.Glink(i).n1 = glink_nodes(2);
                    obj.Glink(i).n2 = glink_nodes(1);
                    % Swap round path order
                    obj.Glink(i).point = fliplr(obj.Glink(i).point);
                end
            end
        end
        
%% HIGH LEVEL METHODS    
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
                disp(['Traversing: Completed ', num2str(i), ' of ', num2str(total_branches)])
            end
        end
        
        
        function obj = FindFWHMall(obj)
            % analyse all airway segments except the trachea.
            disp('Start computing FWHM boundaries of all airway segments')
            total_branches = length(obj.Glink);
            for i = 1:length(obj.Glink)
                % skip the trachea
                %if i == obj.trachea_path
                % TODO: CHANGE THIS!
                if i == 58
                    continue
                end
                obj = FindAirwayBoundariesFWHM(obj, i);
                disp(['FWHMesl: Completed ', num2str(i), ' of ', num2str(total_branches)])
            end
        end
        
%% TRAVERSING AIRWAYS METHODS %%%
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
            
            % The input is the list of ordered index
            % The output is the smooth spline as a matlab sturct
            
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
        
%% TAPERING MEASUREMENT METHODS
        function obj = FindAirwayBoundariesFWHM(obj, link_index)
            %Based on function by Kin Quan 2018 that is based on Kiraly06
            
            slices_size = size(obj.TraversedImage{link_index, 1}, 3);
            
            % Prepping the outputs
            raycast_FWHMl = cell(slices_size, 1);
            raycast_FWHMp = cell(slices_size, 1);
            raycast_FWHMr = cell(slices_size, 1);
            
            % For every traversed slice
            for k = 1:slices_size
                try
                % * Compute airway centre
                % Check that airway centre is slice centre
                center = ...
                    Return_centre_pt_image(...
                    obj.TraversedSeg{link_index, 1}(:,:,k));
                
                % Recompute new centre if necessary
                [centre_ind , new_centre] =  ...
                    Check_centre_with_segmentation(...
                    obj.TraversedSeg{link_index, 1}(:,:,k),...
                    obj.TraversedImage{link_index, 1}(:,:,k));
                if ~centre_ind
                    center = fliplr(new_centre);
                end
                
                % * Raycast
                [CT_rays, seg_rays, coords]= Raycast(obj, ...
                    obj.TraversedImage{link_index, 1}(:,:,k), ...
                    obj.TraversedSeg{link_index, 1}(:,:,k), center);
                
                % * Compute FWHM
                [FWHMl, FWHMp, FWHMr] = AirwaySkel.computeFWHM(CT_rays,...
                    seg_rays, coords);
                
                % * Compute Ellipses
                FWHMl_ellipse = ComputeEllipses(obj, FWHMl);
                FWHMp_ellipse = ComputeEllipses(obj, FWHMp);
                FWHMr_ellipse = ComputeEllipses(obj, FWHMr);
                
                % * Record, catch incase a slice fails.
                    raycast_FWHMl{k,1} = FWHMl_ellipse;
                    raycast_FWHMp{k,1} = FWHMp_ellipse;
                    raycast_FWHMr{k,1} = FWHMr_ellipse;
                catch
                    % TODO: Look at airway 7 when it failed.
                    warning('Fail recorded')
                    raycast_FWHMl{k,1} = NaN;
                    raycast_FWHMp{k,1} = NaN;
                    raycast_FWHMr{k,1} = NaN;                
                end
            end
            obj.FWHMesl{link_index, 1} = raycast_FWHMl;
            obj.FWHMesl{link_index, 2} = raycast_FWHMp;
            obj.FWHMesl{link_index, 3} = raycast_FWHMp;
        end
        
        function [CT_rays, seg_rays, coords] = Raycast(obj, interpslice, interpseg, center)
            % * Compute Rays
            % Getting the range limit of the ray which will be the shortest
            % distance from the centre to the bounadry Need to find the
            % limits of the raw
            image_size = size(interpslice);
            limit_row = abs(image_size(1) - center(1));
            limit_col = abs(image_size(2) - center(2));
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
            
            CT_rays = interp2(interpseg, x_component(:),...
                y_component(:));
            
            seg_rays = interp2(interpslice, x_component(:),...
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
                pi*obj.physical_sampling_interval.^2;
        end
        
        function obj = ComputeTaperValues(obj, link_index)
            obj.specs(link_index).FWHMl_logtaper = ...
                AirwaySkel.ComputeTaperRate(...
                obj.arclength{link_index, 1}, ...
                [obj.FWHMesl{link_index, 1}.area]);
            obj.specs(link_index).FWHMp_logtaper = ...
                AirwaySkel.ComputeTaperRate(...
                obj.arclength{link_index, 1}, ...
                [obj.FWHMesl{link_index, 2}.area]);
            obj.specs(link_index).FWHMr_logtaper = ...
                AirwaySkel.ComputeTaperRate(...
                obj.arclength{link_index, 1}, ...
                [obj.FWHMesl{link_index, 3}.area]);
        end

        function [logtaperrate, cum_arclength, cum_area, path] = ConstructTaperPath(obj, terminal_link_idx)
            G = graph(obj.Gadj);
            % TODO: remove trachea node?
            path = shortestpath(G, obj.carina_node, terminal_link_idx);
            cum_arclength = 0;
            cum_area = [];
            for i = path
                % skip the first of each link except the first one
                if i == path(1)
                    cum_area = [cum_area; obj.FWHMesl{i, 1}{1, 1}.area];
                end
                max_current_arclength = max(cum_arclength);
                current_arclength_array = max_current_arclength+obj.arclength{i,1}(2:end);
                cum_arclength = [cum_arclength; current_arclength_array];

                for j = 2:length(obj.arclength{i,1})
                    try 
                        cum_area = [cum_area; obj.FWHMesl{i, 1}{j, 1}.area];
                    catch
                    	cum_area = [cum_area; NaN];
                    end
                end
            end
            logtaperrate = AirwaySkel.ComputeTaperRate(cum_arclength, cum_area);
            end
            
        
%% VISUALISATION METHODS
        function PlotTree(obj)
            % Plot the airway tree with nodes and links
            % Original Function by Ashkan Pakzad on 27th July 2019.
            
            isosurface(bwskel(obj.seg, 'MinBranchLength', obj.branch_threshold));
            hold on
            
            % edges
            ind = zeros(length(obj.Glink), 1);
            for i = 1:length(obj.Glink)
            ind(i) = obj.Glink(i).point(ceil(end/2));
            end
            [Y, X, Z] = ind2sub(size(obj.seg),ind);
            nums_link = string(1:length(obj.Glink));
            %plot3(X,Y,Z, 'b.', 'MarkerFaceColor', 'none');
            text(X+1,Y+1,Z+1, nums_link, 'Color', [0, 0.3, 0])
            
            % nodes
            X_node = [obj.Gnode.comy];
            Y_node = [obj.Gnode.comx];
            Z_node = [obj.Gnode.comz];
            nums_node = string(1:length(obj.Gnode));
            plot3(X_node,Y_node,Z_node, 'r.', 'MarkerSize', 18, 'Color', 'r');
            text(X_node+1,Y_node+1,Z_node+1, nums_node, 'Color', [0.8, 0, 0])
            
            axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            view(80,0)
            
        end
        
        function Airway3D(obj, link_index)
            % Plot resampled airway slices overlayed with FWHMesl ray cast 
            % points and fitted ellipse
            f = figure('Position',  [100, 100, 850, 600]);
            slide = 1;
            Airway2D(obj, link_index, slide)
            numSteps = size(obj.TraversedImage{link_index,1}, 3);

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
                Airway2D(obj, link_index, val);
            end
            
        end
        
        function Airway2D(obj, link_index, slide)
            % display image
            imagesc(obj.TraversedImage{link_index, 1}(:,:,slide))
            hold on
            colormap gray
            
            try % try block incase FWHMesl has not been executed.
            % plot ray cast results

            plot(obj.FWHMesl{link_index, 1}{slide, 1}.x_points, obj.FWHMesl{link_index, 1}{slide, 1}.y_points,'r.')
            plot(obj.FWHMesl{link_index, 2}{slide, 1}.x_points, obj.FWHMesl{link_index, 2}{slide, 1}.y_points,'c.')
            plot(obj.FWHMesl{link_index, 3}{slide, 1}.x_points, obj.FWHMesl{link_index, 3}{slide, 1}.y_points,'y.')

            % plot ellipse fitting
            ellipse(obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(4),...
                obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(1),...
                obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(2),'m');
            
            ellipse(obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(4),...
                obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(1),...
                obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(2),'b');
            
            ellipse(obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(4),...
                obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(1),...
                obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(2),'y');
            %TODO: set third colour more appropiately
            
            catch
                warning('No FWHMesl data, showing slices without elliptical information.')
            end
            
            % display area measurements
            % TODO: is this needed?
            %dim = [.15 .85 .24 .05];
            %a = annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','y');
            a = rectangle('Position',[0,0,133,10],'FaceColor','y','LineWidth',2);
            ax = gca;
            text(ax, 1,5,sprintf('Arc Length = %4.1f mm; Inner area = %4.2f mm^2; Peak area = %4.2f mm^2; Outer area = %4.2f mm^2; %3.0i of %3.0i', ...
                obj.arclength{link_index, 1}(slide), obj.FWHMesl{link_index, 1}{slide, 1}.area, obj.FWHMesl{link_index, 2}{slide, 1}.area ,obj.FWHMesl{link_index, 3}{slide, 1}.area, slide, size(obj.TraversedImage{link_index, 1},3)));
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
        
        
        function logtaperrate = ComputeTaperRate(arclength, area)        
        p_coeff = polyfit(arclength,log(area),1);
        logtaperrate = p_coeff(1);                   
        end
        
        
        function [FWHMl, FWHMp, FWHMr] = computeFWHM(CT_rays, seg_rays, coords)
            %This is to perfrom the ray casting measurents - the input is in a sturct
            % as
            
            % The basis of the code is base on the description - Virtual Bronchoscopy for Quantitative Airway Analysis
            % by A P Kiraly et al.
            
            %Need to find the length of the ray - this need a loop
            number_of_rays = size(CT_rays,2);
            FWHMl_x = [];
            FWHMl_y = [];
            FWHMp_x = [];
            FWHMp_y = [];
            FWHMr_x = [];
            FWHMr_y = [];
            
            %performing a loop
            for ray = 1:number_of_rays
                
                CT_profile = CT_rays(:,ray);
                seg_profile = seg_rays(:,ray);
                
                % * Find the edge of the interpolated segmentation.
                if seg_profile(1) < 0.5
                    continue % skip if seg edge does not exist
                else
                    % Identify the last vaule that is above the 0.5
                    ind_ray = (seg_profile < 0.5);
                    seg_half = find(ind_ray,1);
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
                
                % TODO: add same anomaly detection as before.
                
                % concat points together
                FWHMl_x = cat(1,FWHMl_x,coords(FWHMl, ray, 1));
                FWHMl_y = cat(1,FWHMl_y,coords(FWHMl, ray, 2));
                FWHMp_x = cat(1,FWHMp_x,coords(FWHMp, ray, 1));
                FWHMp_y = cat(1,FWHMp_y,coords(FWHMp, ray, 2));
                FWHMr_x = cat(1,FWHMr_x,coords(FWHMr, ray, 1));
                FWHMr_y = cat(1,FWHMr_y,coords(FWHMr, ray, 2));
            end
            
            FWHMl = cat(2, FWHMl_x, FWHMl_y);
            FWHMp = cat(2, FWHMp_x, FWHMp_y);
            FWHMr = cat(2, FWHMr_x, FWHMr_y);
        end
    end
end