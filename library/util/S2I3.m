function I = S2I3(sz,I1,I2,I3)
    % Convert sub-indicies to linear indices for 3D arrays.
    %
    % Faster implementation of MATLAB's sub2ind for 3D arrays only. inputs
    % I1, I2 and I3 must be the same length.
    % Reference: https://stackoverflow.com/questions/60954976/fast-general-replacements-for-ind2sub-and-sub2ind-matlab
    % see also sub2ind
    %
    % Args:
    %   sz: size of reference array e.g. `size(3Darray)`
    %   I1: list of corresponding 1st sub-indices
    %   I2: list of corresponding 2nd sub-indices
    %   I3: list of corresponding 3rd sub-indices
    %
    % Return:
    %   I: linear indices of input list

I =  I1 + (I2-1)*sz(1) + (I3-1)*sz(1)*sz(2); 
end