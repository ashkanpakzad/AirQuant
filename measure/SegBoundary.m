classdef SegBoundary < SuperMeasure
    properties
        center
    end
    methods
        function obj = SegBoundary(varargin)
            % See superclass constructor method SuperMeasure.
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
            obj@SuperMeasure(varargin{:});
            assert(~isempty(obj.tube.seg), ['Tube segmentation required for ' ...
                'this measurement method to work'])
        end

        function Measure(obj)
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

            slices_sz = length(obj.tube.seg);
            ellipse = cell(slices_sz, 1);
            ellipse_center = cell(slices_sz, 1);
            
            for k = 1:slices_sz
                try
                    
                    
                    seg_slice = obj.tube.seg{k,1};
                    center =  AirwayCenter(seg_slice, seg_slice);

                    % filter to binary if necessary
                    if length(unique(seg_slice))>2
                        seg_slice = (seg_slice>0.5);
                    end

                    % filter all objects except central
                    label_slice = bwlabel(seg_slice);
                    % get label of central voxel
                    central_l = label_slice(center(1),center(2));
                    if central_l == 0
                        error('central pixel not in segmentation foreground.')
                    end
                    seg_slice_filt = (label_slice == central_l);
                    B = bwboundaries(seg_slice_filt);
                    points = fliplr(B{1,1});
                    
                    % add random jitter to points to avoid caveats in
                    % ellipse
                    jit_magnitude = 1e-3;
                    jitter = unifrnd(-jit_magnitude,jit_magnitude,size(points));
                    points = points+jitter;
                    
                    % compute ellipses
                    ellipse{k,1} = AQEllipse(obj.pixsize, points);
                    ellipse_center{k,1} = center;
                catch e
                    % measure error
                    warning(['SegBoundary Failed: slice ',num2str(k),': ',e.identifier])
                    warning(e.message)
                    ellipse{k,1} = nan;
                    ellipse_center{k,1} = nan;
                end
            end
            obj.measures = ellipse';
            obj.center = ellipse_center;
            
        end

        
    end
end
