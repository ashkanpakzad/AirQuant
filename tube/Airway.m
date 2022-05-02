classdef Airway < Tube
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
    
    properties
        istrachea
        carinaend
    end
    
    methods
        function obj = Airway(inputArg1,inputArg2)
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
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = FWHMMeasurement(obj,inputArg)
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
            outputArg = obj.Property1 + inputArg;
        end
        
    end
end

