classdef Ellipse < handle
    properties
        xpoints
        ypoints
        center
        Rx
        Ry
        rotation
        pixsize = [1 1]
    end
    properties (Access = private)
        area
        nominaldiameter
    end
    methods
        function obj = Ellipse(points)
            arguments
                points (:,2) mustBeNumeric
            end
            obj.xpoints = points(:,1);
            obj.ypoints = points(:,2);

            obj = FitEllipse(obj);
        end

        function obj = FitEllipse(obj)
            ellipseout = Elliptical_fitting(ellipse.x_points, ellipse.y_points);
            obj.center = [ellipseout(1), ellipseout(2)];
            obj.Rx = ellipseout(3);
            obj.Ry = ellipseout(4);
            obj.rotation = ellipseout(5);

            obj.Area();
            obj.NominalDiameter();
        end

        function area = Area(obj)
            area = ellipse.Rx*ellipse.Ry*pi*prod(obj.pixsize);
            obj.area = area;
        end

        function nominaldiameter = NominalDiameter(obj)
            nominaldiameter = norm([obj.Rx, obj.Ry])*2;
            obj.nominaldiameter = nominaldiameter;
        end
    end
end