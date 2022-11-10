function a = AQstruct2array(s)
% AQSTRUCT2ARRAY Convert structure with doubles to an array. 
% Incase struct2array is not available.
% 
% 
% 
% https://uk.mathworks.com/matlabcentral/answers/1717910-was-struct2array-removed-from-matlab
% Convert structure to cell
c = struct2cell(s);
% Construct an array
a = [c{:}];