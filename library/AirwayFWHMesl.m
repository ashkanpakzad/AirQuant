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
            ellipse_peak = cell(slices_sz, 1);
            ellipse_right = cell(slices_sz, 1);
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
                    [innerraypoints, peakraypoints, outerraypoints] = compute_airway_fwhm(...
                        source_rays, seg_rays, coords, outlierremoval);

                    % compute ellipses
                    ellipse_inner{k,1} = AQEllipse(obj.pixsize, innerraypoints);
                    ellipse_peak{k,1} = AQEllipse(obj.pixsize, peakraypoints);
                    ellipse_right{k,1} = AQEllipse(obj.pixsize, outerraypoints);
                    ellipse_center{k,1} = center;
                catch
                    % segmentation exceeds interpolated slice therefore no
                    % measurement recorded.
                    ellipse_inner{k,1} = [];
                    ellipse_peak{k,1} = [];
                    ellipse_right{k,1} = [];
                    ellipse_center{k,1} = [];
                end
            end
            obj.measures = cell(3, length(ellipse_inner));
            obj.measures(1,:) = ellipse_inner;
            obj.measures(2,:) = ellipse_peak;
            obj.measures(3,:) = ellipse_right;
            obj.center = ellipse_center;
            
        end

        function area = OutputArea(obj)
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
            area = cell2mat(cellfun(@(c) [c.area], obj.measures, 'UniformOutput', false));
        end

        function diameter = OutputDiameter(obj)
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
            diameter = cell2mat(cellfun(@(c) [c.diameter], obj.measures, 'UniformOutput', false));
        end
        
    end
end
