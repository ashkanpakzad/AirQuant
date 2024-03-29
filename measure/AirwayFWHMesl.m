classdef AirwayFWHMesl < SuperMeasure
    properties
        center
    end
    methods
        function obj = AirwayFWHMesl(varargin)
            % See superclass constructor method SuperMeasure.
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

            if nargin < 4
                outlierremoval = true;
            end

            slices_sz = length(obj.tube.source);
            ellipse_inner = cell(slices_sz, 1);
            ellipse_outer = cell(slices_sz, 1);
            ellipse_center = cell(slices_sz, 1);
            
            for k = 1:slices_sz
                try
                    
                    
                    source_slice = obj.tube.source{k,1}; 
                    seg_slice = obj.tube.seg{k,1};
                    center =  AirwayCenter(source_slice, seg_slice);

                    [source_rays, seg_rays, coords] = ray_cast(...
                        source_slice, seg_slice, center, ...
                        num_rays, ray_interval);

                    % compute FWHM points of ray
                    [innerraypoints, ~, outerraypoints] = compute_airway_fwhm(...
                        source_rays, seg_rays, coords, outlierremoval);

                    % compute ellipses
                    ellipse_inner{k,1} = AQEllipse(obj.pixsize, innerraypoints);
                    ellipse_outer{k,1} = AQEllipse(obj.pixsize, outerraypoints);
                    ellipse_center{k,1} = center;
                catch
                    % segmentation exceeds interpolated slice therefore no
                    % measurement recorded.
                    ellipse_inner{k,1} = nan;
                    ellipse_outer{k,1} = nan;
                    ellipse_center{k,1} = nan;
                end
            end
            obj.measures = cell(2, length(ellipse_inner));
            obj.measures(1,:) = ellipse_inner;
            obj.measures(2,:) = ellipse_outer;
            obj.center = ellipse_center;
            
        end

        
    end
end
