classdef SuperMeasure < handle
    properties
    tube
    pixsize
    measures
    end
    methods
        function obj = SuperMeasure(varargin)
            obj.pixsize = [varargin{1}, varargin{1}];
            obj.tube = varargin{2};
            obj.Measure(varargin{3:end})
        end
        
        function Measure(obj)
            error('Measure method not set for subclass of SuperMeasure')
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
