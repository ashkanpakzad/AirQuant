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


function spanned_points = ...
    Span_plane_points_from_basis( basis_vector_1, basis_vector_2, ...
    coeff_1 , coeff_2 ,phyical_point )
%Gives all the points that the basis span from the given coeff.

%The input is the basis vector in the from of 1 by 3 (they will be
%resphaped) - The coeff is the array (they will be reshape). Note that the
%coeff need to have the same number if points

%The output is the span points in 3 by n points


%% Need to reshape the variables

basis_vector_1 = basis_vector_1(:);
basis_vector_2 = basis_vector_2(:);
coeff_1 = coeff_1(:)';
coeff_2 = coeff_2(:)';

%% Getting the spanning points
spanned_points = zeros(3,length(coeff_1)*2);

%Need a counter
k = 1;

for i = 1:length(coeff_1)
    for j = 1:length(coeff_2)

        spanned_points(:,k) = coeff_1(i)*basis_vector_1 + ...
            coeff_2(j)*basis_vector_2 + phyical_point;
        k = k+1;
    end
end

end
