function out = CropVol(varargin)
    % Crop volume to limits or non-zero extremes if bool image.
    %
    % Crops to boundary with padding of 10 voxels unless this would
    % overextend beyond the boundary of the source image, then it is
    % cropped to the source image boundary.
    % 
    % TODO
    % ----
    % * Make so lims and cropped source can be output together.
    %
    % Args:
    %   source: source image to be cropped, if only input must be `bool`
    %   lims: `OPTIONAL` extreme limits to crop volume
    % Return:
    %   lims: `OPTIONAL` output if only `source` provided at input.
    %   croppedsource: `OPTIONAL` cropped source image if `lims` not given
    %   at input.
    %   
    s = size(varargin{1});

    if length(varargin) == 1 && islogical(varargin{1})
        % must be a logical input, therefore segmentation
        % find extremes of x, y & z of segmentation and pad by 1.
        all_lin_idx = find(varargin{1}(:) == 1);
        [ix, iy, iz] = ind2sub(s,all_lin_idx);
        out= [
            min(ix) - 10, max(ix) + 10;
            min(iy) - 10, max(iy) + 10;
            min(iz) - 10, max(iz) + 10;
            ];
    elseif length(varargin) == 1 && ~islogical(varargin{1})
        error('If only one input var, expecting logical.')
    elseif length(varargin) == 2
        if size(varargin{2})~= [length(size(varargin{1})),2]
            error('expected 2nd var to be number of dims in 1st var by 2.')
        end
        % incase lower lims are < index 1, set to min possible
        lims = varargin{2};
        for ii = 1:3
            if lims(ii,1) < 1
                lims(ii,1) = 1;
            end
        end
        % incase upper lims are > size dim, set to max possible
        j = 1;
        for ii = 1:3
            if lims(ii,2) > s(j)
                lims(ii,2) = s(j);
            end
            j = j+1;
        end
        out = varargin{1}(lims(1,1):lims(1,2), lims(2,1):lims(2,2), ...
            lims(3,1):lims(3,2));
    end

end
