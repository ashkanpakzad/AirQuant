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

function distance_slice = Generate_distance_slice( slice_of_interest )
%Genrates a map that gives the distence from the centre of the input

%The input is the image slice - this will allow us to dertermine the size
%of the image and the location of the image

%The output is the distance map where each of the voxel

%% Genrating the slice

size_vector = size(slice_of_interest);

%Labelling the disatnce map
distance_slice = false(size_vector(1),size_vector(2));
%Getting the cenre point
centre_pt = Return_centre_pt_image(slice_of_interest);
%labelling the centre
distance_slice(centre_pt(1),centre_pt(2)) = 1;

%Finally perfrom the distenace tarnsfrom
distance_slice = bwdist(distance_slice,'euclidean');

end
