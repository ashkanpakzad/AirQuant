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


function [ind_change,updated_centre] = ...
    Check_centre_with_segmentation( raw_image_slice , seg_image_slice )
%Check if the centre overlays with the actual segmentation - if not then
%the new centre will be based at the centre

%The input is the raw slice and the segmentation slice

%The ouptut is the a bianry indictor, ture means the centre of the slice is
%correct the updated center.

%% Getting the current centre

updated_centre = Return_centre_pt_image(raw_image_slice);
binary_seg = seg_image_slice >= 0.5;

%%
if binary_seg(updated_centre(1),updated_centre(2))
    %No change in the centre
    ind_change = true;
    
else
    %There a change in the centre
    ind_change = false;
    
    %Genrate the 2d distance slice
    distance_slice = Generate_distance_slice(raw_image_slice);
    
    %If the points is missed then the closest connected compomenet is
    %taken frmo the distace tranfrom mask
    index_of_airway = find(binary_seg);
    non_zero_distace_vaules = distance_slice(index_of_airway);
    [~,min_of_min_index] = min(non_zero_distace_vaules);
    
    %Getting the seed point
    seed_point = index_of_airway(min_of_min_index);
    %Need to convert to 2d croods
    seed_point = Return_2d_point_from_index(binary_seg,seed_point);
    
    %% Computing the new centre
    single_compoents_index = ...
        Connected_component_region_2d(seed_point,binary_seg);
    
    %%Need to compute the center of gravity
    [centre_row, centre_col] = ...
        Find_centre_via_distance_transfrom(single_compoents_index);
   
    updated_centre = [centre_row, centre_col];
    
    %Need to be rounded
    updated_centre = round(updated_centre);
    
end

end

