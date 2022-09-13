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
        function obj = Airway(varargin)
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
            obj@Tube(varargin{:});
            obj.istrachea = false;
            obj.carinaend = false;
        end
        
        function obj = SetTrachea(obj)
            obj.istrachea = true;
            obj.generation = 0;
        end

        function obj = SetCarinaEnd(obj)
            obj.SetTrachea();
            obj.carinaend = true;
        end

    end
end

