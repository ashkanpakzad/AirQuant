classdef VesselFWHMesl < SuperMeasure
    properties
        center
    end
    methods
        function obj = VesselFWHMesl(varargin)
            % Measure vessel CT patches directly with aid of the
            % segmentation for expected diameter.
            % 
            % 
            % Based on Quantification of pulmonary vessel diameter in 
            % low-dose CT images by Rudyanto et al. 2015. Modified to 
            % use the segmentation as expected vessel diameter. Identifies
            % peaks with minimum width computed by segmentation in the CT
            % profile.
            %
            %
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            obj@SuperMeasure(varargin{:});
            assert(~isempty(obj.tube.seg), ['Tube segmentation required for ' ...
                'this measurement method to work'])
        end

        function Measure(obj, num_rays, ray_interval, outlierremoval)
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

            if nargin < 2
                num_rays = 30;
            end

            if nargin < 2
                ray_interval = 0.2;
            end

            if nargin < 4
                outlierremoval = true;
            end

            slices_sz = length(obj.tube.source);
            ellipse = cell(slices_sz, 1);
            ellipse_center = cell(slices_sz, 1);
            
            for k = 1:slices_sz
                try
                    
                    
                    source_slice = obj.tube.source{k,1}; 
                    seg_slice = obj.tube.seg{k,1};
                    center =  AirwayCenter(source_slice, seg_slice);

                    [source_rays, seg_rays, coords] = ray_cast2(...
                        source_slice, seg_slice, center, ...
                        num_rays, ray_interval);

                    % compute FWHM points of ray
                    [innerraypoints, ~, outerraypoints] = compute_vessel_fwhm(...
                        source_rays, seg_rays, coords, outlierremoval);

                    % cat array points 
                    arraypoints = cat(1,innerraypoints,outerraypoints);

                    % compute ellipses
                    ellipse{k,1} = AQEllipse(obj.pixsize, arraypoints);
                    ellipse_center{k,1} = center;
                catch e
                    % measure error
                    warning(['VesselFWHM Failed: slice ',num2str(k),': ',e.identifier])
                    warning(e.message)                   
                    ellipse{k,1} = nan;
                    ellipse_center{k,1} = nan;
                end
            end
            obj.measures = ellipse';
            obj.center = ellipse_center;
            
        end

        
    end
end
