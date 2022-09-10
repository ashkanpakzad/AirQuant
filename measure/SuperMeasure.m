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
        
        function var = OutputVar(obj,varname)
            % Collate variable that exists as property of cell of objects
            % in `prop::measures`.
            %
            % Args:
            %   varname(char): variable name to collate. e.g. 'diameter'
            %
            % Return:
            %   var(array)
            %
            try
                var = cell2mat(cellfun(@(c) [c.(varname)], obj.measures, 'UniformOutput', false));
            catch % if nan in measures use for loop
                var = zeros(size(obj.measures));
                for ii = 1:size(obj.measures,2)
                    for jj = 1:size(obj.measures,1)
                        if isa(obj.measures{jj,ii},'AQEllipse')
                            var(jj,ii) = obj.measures{jj,ii}.(varname);
                        else
                            var(jj,ii) = NaN;
                        end
                    end
                end
            end
        end

    end
end
