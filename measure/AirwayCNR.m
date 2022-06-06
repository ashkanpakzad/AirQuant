classdef AirwayCNR < SuperMeasure
    methods
        function obj = AirwayCNR(varargin)
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

        function Measure(obj, modelpath)
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
            
            % set up memory
            slices_sz = length(obj.tube.source);
            ellipse_inner = cell(slices_sz, 1);
            ellipse_outer = cell(slices_sz, 1);
            ellipse_center = cell(slices_sz, 1);
            
            for k = 1:slices_sz
                try
                    % get source slice
                    source_slice = obj.tube.source{k,1};                     
                    
                    % run through model

                    
                    % process inner ellipses
                    innerstruct = struct;
                    innerstruct.center = center;
                    innerstruct.Rx = innerRx;
                    innerstruct.Ry = innerRy;
                    innerstruct.rotation = rotation;
                    
                    % process outer ellipses
                    outerstruct = struct;
                    outerstruct.center = center;
                    outerstruct.Rx = outerRx;
                    outerstruct.Ry = outerRy;
                    outerstruct.rotation = rotation;

                    % compute ellipses
                    ellipse_inner{k,1} = AQEllipse(obj.pixsize, innerstruct);
                    ellipse_outer{k,1} = AQEllipse(obj.pixsize, outerstruct);
                catch
                    % segmentation exceeds interpolated slice therefore no
                    % measurement recorded.
                    ellipse_inner{k,1} = nan;
                    ellipse_outer{k,1} = nan;
                end
            end
            obj.measures = cell(3, length(ellipse_inner));
            obj.measures(1,:) = ellipse_inner;
            obj.measures(3,:) = ellipse_outer;
            
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
