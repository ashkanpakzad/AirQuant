function [output, outvoxdim] = ReorientVolume(source, meta)
    % Reorients input volume from niftiinfo/niftiread into forced LPS
    % orientation by applying necessary axes permutations and flips based
    % on the metadata affine.
    %
    % .. todo :
    %   * make independent of loading mechanism
    %
    % Args:
    %   img: volume from `niftiread`
    %   meta (struct): output from `niftiinfo` for same volume,
    %       containing header information
    % Returns:
    %   [output, outvoxdim]
    %
    %   * output: reoriented volume in LPS configuration
    %   * outvoxdim: `OPTIONAL` output voxel dimensions for LPS reorientation.
    %

    % get affine matrix
    aff_raw_RAS_loaded = meta.Transform.T;
    % set values that are basically 0 to 0
    zero_tol = aff_raw_RAS_loaded < 1e-6 & aff_raw_RAS_loaded > -1e-6;
    aff_raw_RAS = aff_raw_RAS_loaded;
    aff_raw_RAS(zero_tol) = 0;
    if any(aff_raw_RAS_loaded ~= aff_raw_RAS)
        warning('Oblique orientation, precision loss within 1e-6')
    end
    % remove origin information
    aff_raw_RAS(4,1:3) = [0,0,0];
    % remove spacing information
    aff_raw_RAS = aff_raw_RAS./abs(aff_raw_RAS);
    aff_raw_RAS(isnan(aff_raw_RAS)) = 0;
    % convert to LPS affine
    aff_raw_LPS = aff_raw_RAS*diag([-1,-1,1,1]);

    % take absolute values to figure out anatomical axes
    aff_raw_LPS_pos = abs(aff_raw_LPS);

    % identify permutation to achieve L/R,P/A,S/I
    aff_raw_LPS_pos_red = aff_raw_LPS_pos(1:3,1:3);

    newaxes = [1,2,3]; % default
    for i = 1:3
        vec = aff_raw_LPS_pos_red(:,i);
        newaxes(i) = find(vec);
    end

    % identify flips
    aff_raw_LPS_red = aff_raw_LPS(1:3,1:3);
    aff_LPS_red = zeros(3,3);

    % construct new affine with permutation
    for i = 1:3
        aff_LPS_red(:,i) = aff_raw_LPS_red(:,newaxes(i));
    end

    % identfy axes that need flipping
    flips = [0,0,0]; % default
    for i = 1:3
        if sum(aff_LPS_red(:,i)) == -1
            flips(i) = 1;
        end
    end

    % apply permutations and flips
    output = permute(source, newaxes);
    for i = 1:3
        if flips(i) == 1
            output = flip(output, i);
        end
    end

    % Change the meta if requested
    if (nargout > 1)
        v_sz = meta.PixelDimensions;
        outvoxdim = [v_sz(newaxes(1)), v_sz(newaxes(2)),...
            v_sz(newaxes(3))];
    end

end
