classdef SuperMeasure < handle
    properties
    tube
    pixsize
    measures
    end
    methods
        function obj = SuperMeasure(varargin)
            if ~isempty(varargin)
                obj.pixsize = [varargin{1}, varargin{1}];
                obj.tube = varargin{2};
                obj.Measure(varargin{3:end})
            end
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
            try
                area = cell2mat(cellfun(@(c) [c.area], obj.measures, 'UniformOutput', false));
            catch % if nan in measures use for loop
                area = zeros(size(obj.measures));
                for ii = 1:size(obj.measures,2)
                    for jj = 1:size(obj.measures,1)
                        if isa(obj.measures{jj,ii},'AQEllipse')
                            area(jj,ii) = obj.measures{jj,ii}.area;
                        else
                            area(jj,ii) = NaN;
                        end
                    end
                end
            end
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
            try
                diameter = cell2mat(cellfun(@(c) [c.diameter], obj.measures, 'UniformOutput', false));
            catch % if nan in measures use for loop
                diameter = zeros(size(obj.measures));
                for ii = 1:size(obj.measures,2)
                    for jj = 1:size(obj.measures,1)
                        if isa(obj.measures{jj,ii},'AQEllipse')
                            diameter(jj,ii) = obj.measures{jj,ii}.diameter;
                        else
                            diameter(jj,ii) = NaN;
                        end
                    end
                end
            end
        end
        
    end
end
