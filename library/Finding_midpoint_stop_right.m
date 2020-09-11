function index_mid_point = ...
    Finding_midpoint_stop_right( global_ray_profile,...
    mid_threshold , min_idx, max_idx)
%Getting the the mid point of the threshold - assumes the ray profile is
%montonic

% By Ashkan Pakzad, based on function by Kin Quan

%I - global ray profile = CT ray profile
%  - expected FWHM intensity
%  - index of minima
%  - index of maxima
%O - the index of the FWHM

%% Setting up

% get CT intensitys between max and min
index_domain = max_idx:min_idx;
profile_window = global_ray_profile(index_domain);

%Assumes that the inteval is montonic
idx_leq_thresh = mid_threshold <= profile_window;
index_limit = find(idx_leq_thresh,1,'last');

%Now need to find the midpoint - need to consider two cases

if index_limit == 1
    % designated as the maxima.
    index_mid_point = index_domain(1);
    
else
    % multiple corresponding points
    upper_range_int = profile_window(index_limit);
    lower_range_int = profile_window(index_limit + 1);
    
    %Need need to find the closest neigbour
    upper_neig_dist = abs(upper_range_int - mid_threshold);
    lower_neig_dist = abs(lower_range_int - mid_threshold);
    
    if upper_neig_dist <= lower_neig_dist
        index_mid_point = index_domain(index_limit);
    else
        index_mid_point = index_domain(index_limit + 1);
    end
    
end

end

