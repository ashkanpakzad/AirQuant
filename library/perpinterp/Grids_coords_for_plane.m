% Adapted by Ashkan Pakzad (ashkanpakzad.github.io)
%
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

function list_of_coords = ...
    Grids_coords_for_plane( basis_1 , basis_2 , centre_point ,...
    length_of_plane , sampling_interval )
%Genrates that coorindates that we are requires to genreate the intensity

%The input are the 2 basis that span the plane, the sampling interval
%the size of the plane and the orgin point

%The Output will be stuct that contanins the x y z corrds of the plane

%% Getting interpolation points
assert(length_of_plane > 0, 'cannot have negative or 0 plane size')
half_length = length_of_plane/2 - sampling_interval/2;
ceoff_array = single(-half_length:sampling_interval:half_length);

%Getting the basis
list_of_coords = ...
    Span_plane_points_from_basis(basis_1,basis_2,ceoff_array, ceoff_array,...
    single(centre_point));

end

