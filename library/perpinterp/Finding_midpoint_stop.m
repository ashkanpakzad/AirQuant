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

function index_mid_point = ...
    Finding_midpoint_stop( global_ray_profile,...
    mid_threshold , min_limit, max_limit)
%Getting the the mid point of the threshold - assumes the ray profile is
%montonic
%O - the index of the array

%% Setting up

%Assumes that min point comes first
index_domain = min_limit:max_limit;
index_intensity = global_ray_profile(index_domain);

%Assumes that the inteval is montonic
ind_array = mid_threshold <= index_intensity;
index_limit = find(ind_array,1);

%Now need to find the threhold - need to consider two cases
if index_limit == 1

    index_mid_point = index_domain(1);

else

    upper_range_int = index_intensity(index_limit);
    lower_range_int = index_intensity(index_limit - 1);

    %Need need to find the closest neigbour
    upper_neig_dist = abs(upper_range_int - mid_threshold);
    lower_neig_dist = abs(lower_range_int - mid_threshold);

    if upper_neig_dist <= lower_neig_dist
        index_mid_point = index_domain(index_limit);
    else
        index_mid_point = index_domain(index_limit - 1);
    end

end

end
