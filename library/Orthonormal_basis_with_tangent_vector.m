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


function array_of_basis = ...
    Orthonormal_basis_with_tangent_vector( tangent_vector )
%Getting the basis of the vector for the tangent - note that the basis will
%always be in first basis
%Source: https://www.amazon.co.uk/Fundamentals-Computer-Graphics-Peter-Shirley/dp/1568814690
%The input is the tangent vetor - it always must have three enetries

%The ouputs is the basis where the 2nd and 3rd basis is the spans the
%normal vector

%Note that this function has gone through serval modifcations

%% Reshaping the vector

tangent_vector = reshape(tangent_vector,[3,1]);

%Need to conisder the number of zeros
number_of_zeros = sum((tangent_vector == 0));
%Getting the output
array_of_basis = zeros(3,3);

if number_of_zeros < 2
    
    array_of_basis(:,1) = tangent_vector/norm(tangent_vector,2);
    
    %We consider the artibary vector
    array_of_basis(:,2) = [array_of_basis(2,1) -array_of_basis(3,1) 0];
    array_of_basis(:,2) = array_of_basis(:,2)/norm(array_of_basis(:,2),2);
    
    %Getting the last vector
    array_of_basis(:,3) = cross(array_of_basis(:,1),array_of_basis(:,2));
    array_of_basis(:,3) = array_of_basis(:,3)/norm(array_of_basis(:,3),2);
    
    %Perfrom the cross product the the final
    array_of_basis(:,2) = cross(array_of_basis(:,1),array_of_basis(:,3));
    array_of_basis(:,2) = array_of_basis(:,2)/norm(array_of_basis(:,2),2);
else
    
    %Need to find the ith poistion
    nonzero_position = find(tangent_vector);
    
    remaning_position = setdiff([1 2 3] ,nonzero_position);
    
    %Placing the poisition
    array_of_basis(nonzero_position,1) = 1;
    array_of_basis(remaning_position(1),2) = 1;
    array_of_basis(remaning_position(2),3) = 1;
    
end

end

