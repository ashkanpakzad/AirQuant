function [FWHMl, FWHMp, FWHMr] = compute_airway_fwhm(source_rays, seg_rays, ...
        coords, outlier_removal)
    % This is to perfrom the ray casting measurents - the input is in a sturct
    % as
    %
    % The basis of the code is base on the description - Virtual Bronchoscopy for Quantitative Airway Analysis
    % by A P Kiraly et al.

    %Need to find the length of the ray - this need a loop
    number_of_rays = size(source_rays,2);

    FWHMl_r = nan(number_of_rays, 1);
    FWHMp_r = nan(number_of_rays, 1);
    FWHMr_r = nan(number_of_rays, 1);

    %performing a loop
    for ray = 1:number_of_rays

        CT_profile = source_rays(:,ray);
        seg_profile = seg_rays(:,ray);

        % * Find the edge of the interpolated segmentation.
        if seg_profile(1) < 0.5
            continue % skip if seg edge does not exist
        end

        % Identify the last vaule that is above the 0.5
        ind_ray = (seg_profile < 0.5);
        seg_half = find(ind_ray,1);

        if all(ind_ray ~= 1)
            continue % skip if fail to find point <0.5
        end

        % * Find FWHM peak
        [max_int_array , max_location_array] = ...
            findpeaks(CT_profile);

        if isempty(max_int_array)
            continue % skip ray if peak does not exist
        end

        % identify the peak closest to seg half point
        index_diff = abs(max_location_array - seg_half);
        [~,closest_max] = min(index_diff);
        % ensure the point is unique
        FWHMp = max_location_array(closest_max(1));

        % * Compute FWHM on left side
        % loop through profile from FWHM peak to centre
        % find first minima from right to left.
        for i = FWHMp:-1:2
            if ~(CT_profile(i) >= CT_profile(i - 1))
                break
            end
            i = 1; % incase inner is at beginning of profile
        end
        FWHMi = i; % inner FWHM curve

        threshold_int_left = ...
            (CT_profile(FWHMp) + CT_profile(FWHMi))/2;
        FWHMl = Finding_midpoint_stop(CT_profile,threshold_int_left,FWHMi,FWHMp);

        % * Compute FWHM on right side
        % loop through profile from FWHM peak to distal
        % find first minima from left to right.
        for i = FWHMp:length(CT_profile)-1
            if CT_profile(i) <= CT_profile(i + 1)
                break
            end
            i = length(CT_profile);
        end
        FWHMo = i;  % outer FWHM curve
        % TODO: move threshold calculation into function.
        threshold_int_right = ...
            (CT_profile(FWHMp) + CT_profile(FWHMo))/2;
        FWHMr = Finding_midpoint_stop_right(CT_profile,threshold_int_right,FWHMo,FWHMp);

        % concat points together
        FWHMl_r(ray) = FWHMl;
        FWHMp_r(ray) = FWHMp;
        if isempty(FWHMr) % in case right peak not within profile
            FWHMr = NaN;
        end
        FWHMr_r(ray) = FWHMr;

    end

    % anomaly correction on inner ONLY
    valid_rays = [1:number_of_rays]';
    if outlier_removal == true
        [~,anomalies] = rmoutliers(FWHMl_r, 'median');
        valid_rays = valid_rays(~anomalies);
        FWHMl_r = FWHMl_r(~anomalies);
        FWHMp_r = FWHMp_r(~anomalies);
        FWHMr_r = FWHMr_r(~anomalies);
    end

    % convert radial index into plane coordinates
    FWHMl_x = DualIndex(coords(:,:,1), FWHMl_r, valid_rays);
    FWHMl_y = DualIndex(coords(:,:,2), FWHMl_r, valid_rays);
    FWHMp_x = DualIndex(coords(:,:,1), FWHMp_r, valid_rays);
    FWHMp_y = DualIndex(coords(:,:,2), FWHMp_r, valid_rays);
    FWHMr_x = DualIndex(coords(:,:,1), FWHMr_r, valid_rays);
    FWHMr_y = DualIndex(coords(:,:,2), FWHMr_r, valid_rays);

    % output
    FWHMl = cat(2, FWHMl_x, FWHMl_y);
    FWHMp = cat(2, FWHMp_x, FWHMp_y);
    FWHMr = cat(2, FWHMr_x, FWHMr_y);
end