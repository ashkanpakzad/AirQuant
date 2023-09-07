function [FWHM] = compute_vessel_fwhm(source_rays, ...
        seg_rays, coords, outlier_removal)
    % Identify fwhm of vessel using findpeaks on CT image only.

    %Need to find the length of the ray - this need a loop
    number_of_rays = size(source_rays,2);

    FWHMl_r = nan(number_of_rays, 1);
    FWHMp_r = nan(number_of_rays, 1);
    FWHMr_r = nan(number_of_rays, 1);

    %performing a loop
    for ray = 1:number_of_rays
        
        % seg profiling
        S_profile = seg_rays(:,ray);

        % identify expected width based on seg
        midpoint = round(length(S_profile)/2);
        [~,lidx] = min(abs(S_profile(1:midpoint) - 0.5));
        [~,ridx] = min(abs(S_profile(midpoint:end) - 0.5));
        ridx = ridx + midpoint;
        segwidth = abs(ridx-lidx);

        % CT profiling
        CT_profile = source_rays(:,ray);

        % identify peaks, min peak width given by segmentation
        [max_int_array , max_location_array,peakwidth,peakprom] = ...
            findpeaks(CT_profile, MinPeakDistance=segwidth);

        if isempty(max_int_array)
            continue % skip ray if peak does not exist
        end

        % identify the peak closest to the centre
        centre = floor(length(CT_profile)/2);
        index_diff = abs(max_location_array - centre);
        [~,closest_max] = min(index_diff);
        % ensure the point is unique
        FWHMp = max_location_array(closest_max(1));
    
        % find left and right side points at half maximum prominence.
        % peakmax - prom/2
        halfpoint = max_int_array(closest_max(1))-peakprom(closest_max(1))/2;
        belowprom=CT_profile<halfpoint;

        belowpromleft = belowprom(1:FWHMp);
        FWHMl = find(belowpromleft,1,"last");

        belowpromright = belowprom(FWHMp:end);
        FWHMr = find(belowpromright,1,"first")+FWHMp;

        % concat points together
        FWHMl_r(ray) = FWHMl;
        FWHMp_r(ray) = FWHMp;
        FWHMr_r(ray) = FWHMr;

    end

    % anomaly correction
    % concatenate the two sets of right and left
    FWHM_r = cat(1,FWHMl_r,FWHMr_r);
    valid_rays = repmat([1:number_of_rays]',2,1);

    % compute midpoint of profile (center of vessel)
    midpoint = size(source_rays,1)/2;

    if outlier_removal == true
        % identify any outliers by width.
        [~,anomalies] = rmoutliers(abs(FWHM_r-midpoint), 'median');
        valid_rays = valid_rays(~anomalies);
        FWHM_r = FWHM_r(~anomalies);
    end

    % convert radial index into plane coordinates
    FWHM_x = DualIndex(coords(:,:,1), FWHM_r, valid_rays);
    FWHM_y = DualIndex(coords(:,:,2), FWHM_r, valid_rays);

    % output
    FWHM = cat(2, FWHM_x, FWHM_y);
end