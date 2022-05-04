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
        istrachea = false
        carinaend = false
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
        end
        
    end
end

