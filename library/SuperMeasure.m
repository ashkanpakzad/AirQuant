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
            error('Measure method not set for subclass')
        end
        
        function OutputArea(obj)
            error('OutputArea method not set for subclass')
        end

        function OutputDiameter(obj)
            error('OutputDiameter method not set for subclass')
        end
        
    end
end
