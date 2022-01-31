% MIT License
%
% Copyright (c) 2019 Kin Quan
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% By Ashkan Pakzad, ashkanpakzad.github.io based on function by Kin Quan.

function index_mid_point = ...
    Finding_midpoint_stop_right( global_ray_profile,...
    mid_threshold , min_idx, max_idx)
%Getting the the mid point of the threshold - assumes the ray profile is
%montonic


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
