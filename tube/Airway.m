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

        function measureFWHM(obj,outlierremoval)
            % short desc
            %
            % Based on function by Kin Quan 2018 that is based on Kiraly06
            %
            % .. todo: add documentation to this function
            %
            % Args:
            %   x(type):
            %
            % Return:
            %   y(type):
            %
            if nargin < 2
                outlierremoval = true;
            end
            
            slices_sz = length(obj.source);
            ellipse_inner = cell(slices_sz, 1);
            ellipse_peak = cell(slices_sz, 1);
            ellipse_right = cell(slices_sz, 1);
            ellipse_center = cell(slices_sz, 1);
            
            for k = 1:slices_sz
                try
                    center =  AirwayCenter(obj.source{k,1}, obj.seg{k,1});

                    [source_rays, seg_rays, coords] = ray_cast(...
                        obj.source, obj.seg, center, ...
                        obj.num_rays, obj.ray_interval);

                    % compute FWHM points of ray
                    [innerraypoints, peakraypoints, outerraypoints] = compute_airway_fwhm(...
                        source_rays, seg_rays, coords, outlierremoval);

                    % compute ellipses
                    ellipse_inner{k,1} = Ellipse(innerraypoints);
                    ellipse_peak{k,1} = Ellipse(peakraypoints);
                    ellipse_right{k,1} = Ellipse(outerraypoints);
                    ellipse_center{k,1} = center;
                catch
                    % segmentation exceeds interpolated slice therefore no
                    % measurement recorded.
                    ellipse_inner{k,1} = [];
                    ellipse_peak{k,1} = [];
                    ellipse_right{k,1} = [];
                    ellipse_center{k,1} = [];
                end
            end
            obj.patchprop.inner.FWHM.ellipse = ellipse_inner;
            obj.patchprop.peak.FWHM.ellipse = ellipse_peak;
            obj.patchprop.outer.FWHM.ellipse = ellipse_right;
            obj.patchprop.center = ellipse_center;

            CollateMeasures(obj, 'patchprop.inner.FWHM.ellipse.area');
            CollateMeasures(obj, 'patchprop.peak.FWHM.ellipse.area');
            CollateMeasures(obj, 'patchprop.outer.FWHM.ellipse.area');
            CollateMeasures(obj, 'patchprop.inner.FWHM.ellipse.nominaldiameter');
            CollateMeasures(obj, 'patchprop.peak.FWHM.ellipse.nominaldiameter');
            CollateMeasures(obj, 'patchprop.outer.FWHM.ellipse.nominaldiameter');
            
        end
        
        function array = CollateMeasures(obj, propertypattern)
            % collate measures that exist as properties of an object into
            % one array for simplified processing.
            array = [obj.(propertypattern)(:)];
            propsave = strip(propertypattern,'.');
            obj.(propsave) = array;
        end

    end
end

