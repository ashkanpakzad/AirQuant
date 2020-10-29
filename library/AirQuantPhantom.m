classdef AirQuantPhantom < handle % handle class
    properties
        CT % CT image
        CTinfo % CT metadata
        seg % binary airway segmentation
        % CT Properties/resampling params
        max_plane_sz = 40;% max interpolated slice size
        plane_sample_sz = 0.3;% interpolated slice pixel size
        spline_sample_sz = 0.25;% mm interval to sample along branch arclength
        plane_scaling_sz = 5; % scale airway diameter approx measurement.
        % Ray params
        num_rays = 50; % check methods
        ray_interval = 0.2; % check methods
        skel % skeleton based on segementation
        savename % filename to load/save
    end
    properties (SetAccess = private)
        % Graph Properties
        Gadj % undirected Graph Adjacency matrix
        Gnode % Graph node info
        Glink % Graph edge info
        Gdigraph % digraph object
        skel_points % list of skeleton points
        % Resampled image slices along graph paths/airway segments.
        Dmap % distance transform of seg
        Splines % pp splines of branches
        TraversedImage % perpendicular CT slices of all airways
        TraversedSeg % perpendicular segmentation slices of all airways
        arclength % corresponding arclength measurement traversed slices
        FWHMesl % FWHMesl algorithm results for every airway
        Specs % Airway specs
    end
    
    methods
        %% INITIALISATION METHODS
        function obj = AirQuantPhantom(CTimage, CTinfo, segimage, skel, savename, params)
            % Initialise the AirQuant class object.
            % if using default settings, do not provide params argument.
            
            % Can provide just file name, if analyses previously complete.
            if nargin == 1
                savename = CTimage;
                obj.savename = savename;
            end
            
            if isfile(savename)
                load(savename)
                obj.savename = savename;
            else
                obj.CTinfo = CTinfo;
                obj.CT = reorientvolume(CTimage, obj.CTinfo);
                % TODO: consider preprocess segmentation to keep largest
                % connected component.
                % ensure no holes in segmentation
                obj.seg = reorientvolume(imfill(segimage,'holes'), obj.CTinfo);
                % set params
                if nargin > 5
                    obj.max_plane_sz = params.max_plane_sz;
                    obj.plane_sample_sz = params.plane_sample_sz;
                    obj.spline_sample_sz = params.spline_sample_sz;
                    obj.num_rays = params.num_rays;
                    obj.ray_interval = params.ray_interval;
                end
                % graph airway skeleton
                obj = GenerateSkel(obj,skel);
                % Convert into digraph
                obj = AirwayDigraph(obj);
                % Compute distance transform
                obj.Dmap = bwdist(~obj.seg);
                % set up empty cell for traversed CT and segmentation
                obj.Splines = cell(length(obj.Glink),1);
                obj.TraversedImage = cell(length(obj.Glink),1);
                obj.TraversedSeg = cell(length(obj.Glink),1);
                % set up empty specs doubles
                obj.arclength = cell(length(obj.Glink),1);
                obj.Specs = struct();
                % set up empty cell for recording raycast/fwhmesl method
                obj.FWHMesl = cell(length(obj.Glink),3);
                % save class
                obj.savename = savename;
                save(obj)
            end
        end
        
        function obj = GenerateSkel(obj, skel)
            % Generate the airway skeleton
            % if skeleton not provided, generate one.
            if isempty(skel)
                obj.skel = Skeleton3D(obj.seg);
            else
                obj.skel = reorientvolume(skel, obj.CTinfo);
            end
            % create graph from skeleton.
            [obj.Gadj,obj.Gnode,obj.Glink] = Skel2Graph3D(obj.skel,0);
            
        end
        
        function obj = AirwayDigraph(obj)
            % Converts the output from skel2graph into a digraph tree
            % network, such that there are no loops and edges point from
            % the trachea distally outwards.
            
            % Create digraph with edges in both directions, loop through
            % and remove if found later in the BF search.
            G = digraph(obj.Gadj);
            % half of the edges will be removed
            removal = zeros(height(G.Edges)/2 , 1);
            j = 1;
            for i = 1:height(G.Edges)
                % identify the z slice of endnodes
                rank1 = obj.Gnode(G.Edges.EndNodes(i,1)).comz;
                rank2 = obj.Gnode(G.Edges.EndNodes(i,2)).comz;
                if rank2 > rank1
                    % record if not central-->distal
                    removal(j) = i;
                    j = j + 1;
                end
            end
            removal(removal == 0) = [];
            G = rmedge(G,removal);
            
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
            
            % Copy over Glink to create digraph with same edge order
            edges = [[obj.Glink.n1]', [obj.Glink.n2]'];
            weights = zeros(length(obj.Glink),1);
            for i = 1:length(obj.Glink)
                weights(i) = length(obj.Glink(i).point);
            end
            labels = [1:length(obj.Glink)]';
            Edgetable = table(edges,weights,labels,'VariableNames',{'EndNodes', 'Weight', 'Label'});
            
            obj.Gdigraph = digraph(Edgetable);
            % add node properties from Gnode
            obj.Gdigraph.Nodes.comx(:) = [obj.Gnode(:).comx];
            obj.Gdigraph.Nodes.comy(:) = [obj.Gnode(:).comy];
            obj.Gdigraph.Nodes.comz(:) = [obj.Gnode(:).comz];
            obj.Gdigraph.Nodes.ep(:) = [obj.Gnode(:).ep];
            obj.Gdigraph.Nodes.label(:) = [1:length(obj.Gnode)]';
            
        end
        
        %% UTILITIES
        
        function save(obj)
            save(obj.savename, 'obj', '-v7.3')
        end
        
        function obj = ComputeSkelPoints(obj)
            [XP, YP, ZP] = ind2sub(size(obj.seg), find(obj.skel == 1));
            obj.skel_points = [XP, YP, ZP]; % list of skel points
        end
        
        function branch_seg = ClassifySegmentation(obj)
            % TODO: consider making this function more robust!
            % label every branch in segmentation to AS branch index.
            
            % Find linear indicies of skeleton
            skel_ind = find(obj.skel == 1);
            classed_skel = zeros(size(skel_ind));
            % for each skeleton branch
            for j = 1:length(obj.Glink)
                % get list of points in branch
                for m = 1:length(obj.Glink(j).point)
                    % put edge number by skeleton index
                    ind = find(skel_ind == obj.Glink(j).point(m));
                    classed_skel(ind) = j;
                end
            end
            
            % get list of everypoint on segmentation
            [XPQ, YPQ, ZPQ] = ind2sub(size(obj.seg),find(obj.seg == 1));
            PQ = [XPQ,YPQ,ZPQ];
            % get list of everypoint on skeleton
            if isempty(obj.skel_points)
                ComputeSkelPoints(obj)
            end
            P = obj.skel_points;
            % find nearest seg point to it on skeleton
            T = delaunayn(P);
            k = dsearchn(P,T,PQ);
            % find that skeleton point's edge assignment
            branch_seg = zeros(size(obj.seg));
            for i = 1:length(PQ)
                branch_seg(PQ(i,1),PQ(i,2),PQ(i,3)) = classed_skel(k(i));
            end
        end
        
        function report = debuggingreport(obj)
            % produce a table showing success of processing for each airway
            % branch.
            report = struct('airway',num2cell(1:length(obj.Glink)));
           
            % check for each airway arclength and FWHM failures
            for i = 1:length(obj.Glink)
                report(i).arclength = ~any(isnan(obj.arclength{i})) & ~any(isempty(obj.arclength{i}));
                try
                report(i).FWHM_inner = ~isempty(obj.FWHMesl{i,1});
                report(i).FWHM_peak = ~isempty(obj.FWHMesl{i,2});
                report(i).FWHM_outer = ~isempty(obj.FWHMesl{i,3});
                catch
                end
            end
                        % remove trachea branch. 
            report(obj.trachea_path) = [];
            report = struct2table(report);
        end
                  
        function report = characteristicsreport(obj, type, showfig)
            % produce a table showing 
            % number of airways per generation
            % number of airways per generation per lobe
            rawgen = [obj.Glink(:).generation];
            generations = unique(rawgen);
%             gencount = zeros(size(generations));
%             for i = 1:length(generations)
%                 for j = 1:length(obj.Glink)
%                     if obj.Glink(j).generation == generations(i)
%                     end
%                 end
%             end
            [gencount, ~] = histcounts(rawgen,generations);

            switch type
                case 'generation'
                    report = table(generations, gencount);
                case 'lobe'
                    
                otherwise
                    error('Choose type generation or lobe')
            end
            
        end
        %% HIGH LEVEL METHODS
        function obj = AirwayImageAll(obj)
            % Traverse all airway segments except the trachea.
            disp('Start traversing all airway segments')
            total_branches = length(obj.Glink);
            % check to see if any branches already processed.
            incomplete = cellfun(@isempty, obj.TraversedImage);
            
            for i = 1:length(obj.Glink)
                % skip the trachea or already processed branches
                obj = CreateAirwayImage(obj, i);
                disp(['Traversing: Completed ', num2str(i), ' of ', num2str(total_branches)])
            end
            disp('Traversing: Done')
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
                obj = FindAirwayBoundariesFWHM(obj, i);
                disp(['FWHMesl: Completed ', num2str(i), ' of ', num2str(total_branches)])
            end
            save(obj)
            disp('FWHMesl: Done')
        end
        
        %% SPLINE METHODS
        function ComputeSpline(obj, link_index)
            % Computes a smooth spline of a single graph edge.
            % Based on original function by Kin Quan 2018
            
            % The input is the list of ordered index
            % The output is the smooth spline as a matlab sturct
            
            %Convert into 3d Points
            [x_point, y_point, z_point] = ind2sub(size(obj.CT),obj.Glink(link_index).point);
            
            %Smoothing the data
            voxel_sz = obj.CTinfo.PixelDimensions;
            smooth_x = smooth(x_point*voxel_sz(1));
            smooth_y = smooth(y_point*voxel_sz(2));
            smooth_z = smooth(z_point*voxel_sz(3));
            
            %Complete smooth data
            smooth_data_points = [smooth_x smooth_y smooth_z]';
            
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
            
            obj.arclength{link_index, 1} = arc_length;
            
    end
        
        %% TRAVERSING AIRWAYS METHODS %%%
        function obj = CreateAirwayImage(obj, link_index)
            % Constructs perpendicular images as if travelling along an
            % airway segment in CT image and Segmentation.
            spline_t_points = ComputeSplinePoints(obj, link_index);
            % set up slice store
            TransAirwayImage = cell(length(spline_t_points),1);
            TransSegImage = cell(length(spline_t_points),1);
            for i = 1:length(spline_t_points)
%                 try
                    spline = obj.Splines{link_index, 1};
                    % * Compute Normal Vector per spline point
                    [normal, CT_point] = AirQuant.ComputeNormal(spline, ...
                        spline_t_points(i));
                    % get approx airway size from distance map
                    approx_diameter = ComputeDmapD(obj, CT_point);
                    % compute intepolated slice size
                    plane_sz = ceil(approx_diameter*obj.plane_scaling_sz);
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
%                 catch
%                     % TODO identify why arc_length cannot be calculated in
%                     % some cases.
%                     warning('Failed to interpolate slice/identify arclength')
%                     TransAirwayImage(:,:,i) = zeros(133);
%                     TransSegImage(:,:,i) = zeros(133);
%                     obj.arclength{link_index, 1}(i) = NaN;
%                 end
            end
            % * Save traversed image and arclength for each image
            obj.TraversedImage{link_index, 1} = TransAirwayImage;
            obj.TraversedSeg{link_index, 1} = TransSegImage;
            % save obj to disk after every branch
            save(obj)
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
            % TODO: Look at what these two lines do.
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
        function obj = FindAirwayBoundariesFWHM(obj, link_index)
            %Based on function by Kin Quan 2018 that is based on Kiraly06
            
            slices_sz = size(obj.TraversedImage{link_index, 1}, 1);
            
            % Prepping the outputs
            raycast_FWHMl = cell(slices_sz, 1);
            raycast_FWHMp = cell(slices_sz, 1);
            raycast_FWHMr = cell(slices_sz, 1);
            
            % For every traversed slice
            for k = 1:slices_sz
                try
                    % * Compute airway centre
                    % Check that airway centre is slice centre
                    center = ...
                        Return_centre_pt_image(...
                        obj.TraversedSeg{link_index, 1}{k,1});
                    
                    % Recompute new centre if necessary
                    [centre_ind , new_centre] =  ...
                        Check_centre_with_segmentation(...
                        obj.TraversedImage{link_index, 1}{k,1}, ...
                        obj.TraversedSeg{link_index, 1}{k,1});
                    if ~centre_ind
                        center = fliplr(new_centre);
                    end
                    
                    % * Raycast
                    [CT_rays, seg_rays, coords]= Raycast(obj, ...
                        obj.TraversedImage{link_index, 1}{k,1}, ...
                        obj.TraversedSeg{link_index, 1}{k,1}, center);
                    
                    % * Compute FWHM
                    [FWHMl, FWHMp, FWHMr] = AirQuant.computeFWHM(CT_rays,...
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
                    % segmentation exceeds interpolated slice therefore no
                    % measurement recorded.
                    raycast_FWHMl{k,1} = NaN;
                    raycast_FWHMp{k,1} = NaN;
                    raycast_FWHMr{k,1} = NaN;
                end
            end
            obj.FWHMesl{link_index, 1} = raycast_FWHMl;
            obj.FWHMesl{link_index, 2} = raycast_FWHMp;
            obj.FWHMesl{link_index, 3} = raycast_FWHMr;
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
     
        %% TAPERING ANALYSIS METHODS
        
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
            terminalnodelist(terminalnodelist ==  obj.trachea_node) = [];
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
        
        function [intrataper, averagearea] = ComputeIntraTaperAll(obj, prunelength)
            % loop through branches
            intrataper = NaN(length(obj.arclength), 3);
            averagearea = NaN(length(obj.arclength), 3);
            for ii = 1:length(obj.arclength)
                if ii == obj.trachea_path
                    continue
                end
                
                [intrataper(ii,:), averagearea(ii,:)] = ComputeIntraTaper(obj, prunelength, ii);
            end
        end
        
        function [intrataper, averagearea] = ComputeIntraTaper(obj, prunelength, idx, plotflag)
            % prunelength given as a 2 element vector, the length in mm to
            % ignore at the begining and end of the branch.
            
            if nargin < 4
                plotflag = 0;
            end
            
            intrataper = NaN(1, 3);
            averagearea = NaN(1, 3);
            
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
            disp(idx)
            prune = (al >= prunelength(1) & al <= al(end) - prunelength(2));
            al = al(prune);
%             areas = areas(repmat(prune',1,3));
            coeff = NaN(2,3);
            % convert area to radii
            areas = sqrt(areas/pi);
            for jj = 1:3
                try % incase no branch left after pruning/too few points
                    areavec = areas(prune,jj);
                    % fit bisquare method
                    coeff(:,jj) = robustfit(al, areavec,'bisquare');
                    % compute intra-branch tapering as percentage
                    intrataper(jj) = -coeff(2,jj)/coeff(1,jj) * 100;
                    % compute average area
                    averagearea(jj) = mean(areavec, 'omitnan');
                catch
                    % leave as NaN
                end
            end
            
            if plotflag == 1 && ~any(isnan(intrataper))
                titlevec = ["inner"; "peak"; "outer"];
                for jj = 1:3
                subplot(3,1,jj)
                plot(al, areas(prune,jj), 'k.')
                hold on
                plot(al, coeff(2,jj)*al + coeff(1,jj),'r')
                legend('data', 'bisquare fit', 'Location', 'best')
                xlabel('arclength (mm)')
                ylabel('area mm^2')
                title(sprintf('Branch idx: %i %s intrataper value: %.2f%% average: %.2f mm.',...
                idx, titlevec(jj), intrataper(jj), averagearea(jj)))
                hold off
                
                end
            end
            
            end
        
        function intertaper = ComputeInterTaper(obj, averagearea)
            % use output from intrataperall
        % loop through branches
            intertaper = NaN(length(obj.arclength), 3);
            for ii = 1:length(averagearea)
                if ii == obj.trachea_path
                    continue
                end
                for jj = 1:3
                    % identify parent by predecessor node
                    parent = find([obj.Glink.n2] == obj.Glink(ii).n1);
                    intertaper(ii,jj) = (averagearea(parent, jj) - averagearea(ii,jj))...
                        /(averagearea(parent, jj)) * 100;
                end
            end

        end
        
        function SegmentTaperResults = SegmentTaperAll(obj, prunelength)
            
            % compute taper results by segment
            [intrataper, avg] = ComputeIntraTaperAll(obj, prunelength);
            intertaper = ComputeInterTaper(obj, avg);
            
            % organise into column headings
            branch = 1:length(obj.Glink);
            
            inner_intra = intrataper(:, 1);
            peak_intra = intrataper(:, 2);
            outer_intra = intrataper(:, 3);
            
            inner_avg = avg(:, 1);
            peak_avg = avg(:, 2);
            outer_avg = avg(:, 3);
            
            inner_inter = intertaper(:, 1);
            peak_inter = intertaper(:, 2);
            outer_inter = intertaper(:, 3);
            
            % convert to table
            SegmentTaperResults = table(branch', inner_intra, peak_intra, ...
                outer_intra, inner_avg, peak_avg, outer_avg, ...
                inner_inter, peak_inter, outer_inter);
            
            % add gen info
            SegmentTaperResults.generation = [obj.Glink.generation]';
            
            % add lobe info if available
            try % only add lobe information if it exists
                SegmentTaperResults.lobe = [obj.Glink.lobe]';
                % sort by lobe
                SegmentTaperResults = sortrows(SegmentTaperResults, 'lobe');
            catch
            end
            
            % Save to AQ object
            obj.Specs.SegmentTaperResults = SegmentTaperResults;
        end
        
        %% TAPERING VISUALISATION METHODS
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
        
        %% VISUALISATION METHODS
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
%             h = plot(G, 'Layout', 'layered');
            h.NodeColor = 'r';
            h.EdgeColor = 'k';
        end
        
        function h = plotgenlabels(obj)
            gens = [AS.Glink(:).generation];
            genslabels = gens(obj.Gdigraph.Edges.Label);
            h = plot(obj, genslabels);
        end
        
        function PlotTree(obj)
            % Plot the airway tree with nodes and links in image space. Set
            % gen to the maximum number of generations to show.
            % Original Function by Ashkan Pakzad on 27th July 2019.

            
%             isosurface(obj.skel);
%             alpha(0.7)
            
            isosurface(obj.skel)
            alpha(0.7)
            
            hold on
            
            % edges
            ind = zeros(length(obj.Glink), 1);
            for i = 1:length(obj.Glink)
                ind(i) = obj.Glink(i).point(ceil(end/2));
            end
            [Y, X, Z] = ind2sub(size(obj.skel),ind);
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
            
            %axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            view(80,0)
            axis vis3d
            % undo matlab display flip
            ax = gca;
            ax.XDir = 'reverse';
            
        end
        
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
        end
        
        function PlotMap3D(obj, mode)
            % mode = 'TaperGradient', 'generation', 'lobes'
            axis([0 size(obj.CT, 1) 0 size(obj.CT, 2) 0 size(obj.CT, 3)])
            
            % generating the color data
            cdata = zeros(size(obj.seg));
            branch_seg = ClassifySegmentation(obj);
            switch mode 
                case 'TaperGradient'
                    % TODO: rewrite this bit....
                    for i = 1:length(obj.Specs)
                        cdata(branch_seg == i) = obj.Specs(i).FWHMl_logtaper*-1;
                    end
                case 'Generation'
                    for i = 1:length(obj.Glink)
                        cdata(branch_seg == i) = obj.Glink(i).generation;
                    end
                    clims = [0 max(cdata(:))];
%                     colorbarstring = 'Generation Number';
                case 'Lobe'
                    % convert lobe id to number
                    lobeid = {'B','RU','RM','RL','LU','LUlin','LL'};
                    for i = 1:length(obj.Glink)
                        cdata(branch_seg == i) = find(strcmp(lobeid, obj.Glink(i).lobe))-1;
                    end
                    clims = [0 max(cdata(:))+1];
%                     colorbarstring = 'Lobe';
%                     colourshow = clims(1):clims(2);
%                     colourlabels = lobeid;
                    
                otherwise
                    error('Choose appropiate mode.')
            end
            % producing segmentation 3d object
            p = patch(isosurface(obj.seg));
            % map colour indices to 3d object vertices
            isocolors(cdata, p);
            % % Reassign colors to meaningful values
            Vertcolour = p.FaceVertexCData;
            vcolours = uniquetol(Vertcolour);
            colourind = unique(cdata);
            for i = 1:length(vcolours)
                Vertcolour(Vertcolour==vcolours(i)) = colourind(i);
            end
            % overwrite vertices colour with meaningful values
            p.FaceVertexCData = Vertcolour;
            
            p.FaceColor = 'interp';
            p.EdgeColor = 'none';
            map = linspecer(max(cdata(:))+2);
            colormap(map)
%             c = colorbar('Ticks', colourshow, 'TickLabels', colourlabels);
%             c.Label.String = colorbarstring;
            caxis(clims)
            view(80,0)
            axis vis3d
            % undo matlab display flip
            ax = gca;
            ax.XDir = 'reverse';
        end
        
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
            
            try % try block incase FWHMesl has not been executed.
                % plot ray cast results
                
%                 plot(obj.FWHMesl{link_index, 1}{slide, 1}.x_points, obj.FWHMesl{link_index, 1}{slide, 1}.y_points,'r.')
%                 plot(obj.FWHMesl{link_index, 2}{slide, 1}.x_points, obj.FWHMesl{link_index, 2}{slide, 1}.y_points,'c.')
%                 plot(obj.FWHMesl{link_index, 3}{slide, 1}.x_points, obj.FWHMesl{link_index, 3}{slide, 1}.y_points,'y.')
                
                % plot ellipse fitting
                ellipse(obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(4),...
                    obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(1)+min_centre,...
                    obj.FWHMesl{link_index, 1}{slide, 1}.elliptical_info(2)+min_centre,'m');
                
%                 ellipse(obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(4),...
%                     obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(1),...
%                     obj.FWHMesl{link_index, 2}{slide, 1}.elliptical_info(2),'b');
                
                ellipse(obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(3),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(4),...
                    obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(5),obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(1)+min_centre,...
                    obj.FWHMesl{link_index, 3}{slide, 1}.elliptical_info(2)+min_centre,'y');
                %TODO: set third colour more appropiately
                
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
        
        %% EXPORT METHODS
        function exportlobes(obj, savename)
            % export airway segmentation labelled by lobes to nii.gz
            
            % classify by branch first
            branch_seg = ClassifySegmentation(obj);
            lobeid = {'B','RU','RM','RL','LU','LUlin','LL'};
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
            header.Description = 'Airway Lobe Segmentation using AirQuant, layers: B,RU,RM,RL,LU,LUlin,LL';
            
            niftiwrite(lobe_airway_seg, savename, header, 'Compressed', true);
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
%                 elseif seg_profile(length(seg_profile)) > 0.5
%                     % set to slice edge if segmentation exceeds slice
%                     seg_half = length(ind_ray);
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