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


function arclength_distance = Arc_length_to_point(matlab_spline_sturct, ...
    para_stop)
%Find the arclength at from the start to the given stopping parameter.

%The input is the point of the parameterization and the spline sturct in
%matlab

%The ouput is the arclength

%% Setting up the inputs

%Setting up the points

%Taking the sturcts apart
[break_points , ~, ~ ,~,~] = ppbrk(matlab_spline_sturct);

%Finding the correct break point to stop;
limit_break_point = (break_points > para_stop);
number_of_pieces = find(limit_break_point,1);

%% Getting the ingtral point.

%Need to get the differnial array
diff_struct = fnder(matlab_spline_sturct);
[~,coeff_array,~,~] = ppbrk(diff_struct);

%Settinig up the loop
total_arclength = 0;

for i = 1:(number_of_pieces - 2)

    %We will be intergarting each of the spline speratlate;

    %Consider the first interval
    interval_upper_bound = break_points(i+1) - break_points(i);

    %Getting the coeff
    coeff_index = 3*i - [2 1 0];
    coeff_array_index = coeff_array(coeff_index,:);
    %Setting the integral function
    intergate_function = @(x) sqrt( ...
        (coeff_array_index(1,1)*x.^2 + coeff_array_index(1,2)*x + coeff_array_index(1,3)).^2 + ...
        (coeff_array_index(2,1)*x.^2 + coeff_array_index(2,2)*x + coeff_array_index(2,3)).^2 + ...
        (coeff_array_index(3,1)*x.^2 + coeff_array_index(3,2)*x + coeff_array_index(3,3)).^2 );


    %Perfroming the intergration
    intergated_interval = integral(intergate_function,0,interval_upper_bound);

    total_arclength = total_arclength + intergated_interval;

end

%% Finding the incompetle interval

final_index = number_of_pieces - 1;
final_interval = para_stop - break_points(final_index);

%Getting the coeff
coeff_index = 3*(final_index) - [2 1 0];
coeff_array_index = coeff_array(coeff_index,:);
%Setting the integral function
intergate_function = @(x) sqrt( ...
    (coeff_array_index(1,1)*x.^2 + coeff_array_index(1,2)*x + coeff_array_index(1,3)).^2 + ...
    (coeff_array_index(2,1)*x.^2 + coeff_array_index(2,2)*x + coeff_array_index(2,3)).^2 + ...
    (coeff_array_index(3,1)*x.^2 + coeff_array_index(3,2)*x + coeff_array_index(3,3)).^2 );

arclength_distance = total_arclength + ...
    integral(intergate_function,0,final_interval);




end
