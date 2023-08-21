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
        min_diameter
        maj_diameter
        diameter
        hydraulic_diameter
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
            
            % compute values from ellipses
            if obj.Rx < 0 || obj.Ry < 0 || isnan(obj.Rx) || isnan(obj.Ry)
                % if ellipse RA or RB turns out to be negative
                obj.area = nan;
                obj.hydraulic_diameter = nan;
                obj.diameter = nan;
                obj.min_diameter = nan;
                obj.maj_diameter = nan;
            else
                obj.Area();
                obj.MinMajDiameters();
                obj.NominalDiameter();
                obj.HydraulicDiameter();
            end
        end

        function obj = FitEllipse(obj)
            points = unique([obj.xpoints, obj.ypoints],'rows');
            if size(points) < 5
                error('Need at least 5 points to fit ellipse.')
            end
            ellipseout = Elliptical_fitting(points(:,1), points(:,2));
            obj.center = [ellipseout(1), ellipseout(2)];
            obj.Rx = ellipseout(3);
            obj.Ry = ellipseout(4);
            obj.rotation = ellipseout(5);
        end

        function area = Area(obj)
            area = obj.Rx*obj.Ry*pi*prod(obj.pixsize);
            obj.area = area;
        end

        function [min_diameter,maj_diameter] = MinMajDiameters(obj)
            radiusx = obj.Rx*obj.pixsize(1);
            radiusy = obj.Ry*obj.pixsize(1);
            min_diameter = min(radiusx,radiusy)*2;
            maj_diameter = max(radiusx,radiusy)*2;
            obj.min_diameter = min_diameter;
            obj.maj_diameter = maj_diameter;
        end

        function nominaldiameter = NominalDiameter(obj)
            nominalradius = sqrt(obj.Rx*obj.Ry)*obj.pixsize(1);
            nominaldiameter = nominalradius*2;
            obj.diameter = nominaldiameter;
        end

        function hydraulic_diameter = HydraulicDiameter(obj)
            % Compute the hydraulic diameter for this ellipse.
            % 
            % The hydraulic diameter is of mechanical interest, computed by
            % 4*area/perimeter of the cross section. It should be noted
            % that the area is unlikely to be a uniform cross section, i.e.
            % not a circle.
            %
            % Perimeter of an ellipse is not a straightforward calculation.
            % This method uses 'ellipsePerimeter' by Santiago Benito. It
            % uses Infinite series expansion, up to the order of 5. 
            %
            %

            % get/compute area
            if ~isempty(obj.area)
                area = obj.area;
            else
                area = obj.Area();
            end

            % compute perimeter
            perimeter = ellipsePerimeter(obj.Rx,obj.Ry,5);

            % compute
            hydraulic_diameter = 4*area/perimeter;
            obj.hydraulic_diameter = hydraulic_diameter;
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