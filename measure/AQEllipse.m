classdef AQEllipse < handle
    properties
        xpoints = []
        ypoints = []
        center
        Rx
        Ry
        rotation
        pixsize
    end
    properties (SetAccess = protected)
        area
        diameter
    end
    methods
        function obj = AQEllipse(pixsize, varargin)
            obj.pixsize = pixsize;
            if nargin == 2 && ~isa(varargin{1},"struct")
                % input is points only
                points = varargin{1};
                assert(any(size(points)==2), 'Expected one argument input to be points of size [2 :] or [: 2].')
                obj.xpoints = points(:,1);
                obj.ypoints = points(:,2);

                obj = FitEllipse(obj);
            elseif nargin == 2 && isa(varargin{1},"struct")
                % input is a struct with reference to ellipse props
                obj.center = varargin{1}.center;
                obj.Rx = varargin{1}.Rx;
                obj.Ry = varargin{1}.Ry;
                obj.rotation = varargin{1}.rotation;

            elseif nargin == 6
                % input is a vector with ellipse props in order
                % [centerx, centery, Rx, Ry, rotation]
                obj.center = [varargin{1}, varargin{2}];
                obj.Rx = varargin{3};
                obj.Ry = varargin{4};
                obj.rotation = varargin{5};
            end
            
            obj.Area();
            obj.NominalDiameter();
        end

        function obj = FitEllipse(obj)
            ellipseout = Elliptical_fitting(obj.xpoints, obj.ypoints);
            obj.center = [ellipseout(1), ellipseout(2)];
            obj.Rx = ellipseout(3);
            obj.Ry = ellipseout(4);
            obj.rotation = ellipseout(5);
        end

        function area = Area(obj)
            area = obj.Rx*obj.Ry*pi*prod(obj.pixsize);
            obj.area = area;
        end

        function nominaldiameter = NominalDiameter(obj)
            nominaldiameter = norm([obj.Rx, obj.Ry])*2;
            obj.diameter = nominaldiameter;
        end

        function [he, hp] = plot(obj, min_centre, ax, showellipse, showpoints)
            if nargin < 3
                ax = gca;
            end

            if nargin < 4
            showellipse = true;
            end

            if nargin < 5
            showpoints = false;
            end
            
            if showellipse == true
                he = ellipse(obj.Rx,obj.Ry,obj.rotation,obj.center(1)+min_centre,...
                        obj.center(2)+min_centre, [], 300, ax);
                he.Color = 'r';
            else
                he = [];
            end

            if showpoints == true
                hold(ax,'on')
                hp = plot(ax, obj.xpoints+min_centre, obj.ypoints+min_centre, '.');
                hold(ax,'off')
            else
                hp = [];
            end

        end
    end
end