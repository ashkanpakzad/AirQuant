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

function coeff_to_span = ...
    Construct_coeff_of_spanning_basis( sampling_interval,length_of_plane )
%Getting the ordered coeffient for the possiable linear combination - i.e.
%the coeffient that span the plane.

%The input is the smapling interval (this is usally the voxel size) and the
%length of the plane

%The output will be the array of coeffients that make the spanning set

%% Getting the limits

%Making sure the length of the plane was positive

length_of_plane = ceil(length_of_plane/2);
length_of_plane = abs(length_of_plane);

%Getting the limits
right_array = 0:sampling_interval:length_of_plane;
left_array = 0:-sampling_interval:-length_of_plane;

%Need to reorentate the left array
left_array = fliplr(left_array);

%Joining them together
coeff_to_span = cat(2,left_array,right_array(2:end));

end
