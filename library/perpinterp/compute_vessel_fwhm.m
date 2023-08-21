function [FWHMl, FWHMp, FWHMr] = compute_vessel_fwhm(source_rays, ...
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

    % anomaly correction on inner ONLY
    valid_rays = [1:number_of_rays]';
    if outlier_removal == true
        % identify any outliers by width.
        [~,anomalies] = rmoutliers(FWHMr_r-FWHMl_r, 'median');
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