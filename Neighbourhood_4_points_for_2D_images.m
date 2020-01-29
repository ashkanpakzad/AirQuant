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

function neighbour_points_list_of = ...
    Neighbourhood_4_points_for_2D_images( two_d_points, binary_slice )
%This gives the neigbouring points on a 2D images - we will assume a four
%neighbouring points

%The inputs are in the ini the form of n by 2 and the binary slice i.e.
%where the image come from

%The output is the list of possiable neig

%% Getting the neigiborinig components

%Desingating the output

neighbour_points_list_of = [];
size_of_binary_slice = size(binary_slice);
size_of_input_array = size(two_d_points);
number_of_seeds = size_of_input_array(1);

%% Genrating the neighbourhood

for displacement_in = [-1,1]

    for component = [1,2]

        %Genrating the displacement vector
        displacement_vector = [0 0];
        displacement_vector(component) = displacement_in;
        %need to repurduce the vector for all the points
        displacement_vector = repmat(displacement_vector,number_of_seeds,1);

        %Genrate the neigbours
        neigbouring_points = two_d_points + displacement_vector;

        %Place in output
        neighbour_points_list_of = cat(1, neighbour_points_list_of,...
            neigbouring_points);

    end
end

%% Applying the idicator and pruning the list for repeats

%Removing the inputs
[neighbour_points_list_of, ~] =...
    setdiff(neighbour_points_list_of,two_d_points,'rows');

%Pruning the inputs and applying the function
neighbour_points_list_of = unique(neighbour_points_list_of,'rows');


%% Applying the indicator

size_of_neighbour = size(neighbour_points_list_of);
number_of_neighbour = size_of_neighbour(1);

%Need to genrate the cooditions to st to remove points that are not
%suitable
upper_bound_condition = repmat(size_of_binary_slice,number_of_neighbour,1);
lower_bound_condition = repmat([0 0],number_of_neighbour,1);

%Getting into the inductors

indictor_of_correct_phyical_coords = ((neighbour_points_list_of <= upper_bound_condition) & ...
    (neighbour_points_list_of > lower_bound_condition));

indictor_of_correct_phyical_coords = indictor_of_correct_phyical_coords(:,1)&...
    indictor_of_correct_phyical_coords(:,2);

index_of_correct_points = find(indictor_of_correct_phyical_coords);

neighbour_points_list_of = neighbour_points_list_of(...
    index_of_correct_points,:);


end
