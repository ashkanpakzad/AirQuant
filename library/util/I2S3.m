function [I1, I2, I3] = I2S3(sz,I)
    % Convert linear indicies to sub-indices for 3D arrays.
    %
    % Faster implementation of MATLAB's ind2sub for 3D arrays only.
    % Reference: https://stackoverflow.com/questions/60954976/fast-general-replacements-for-ind2sub-and-sub2ind-matlab
    % see also ind2sub
    %
    % Args:
    %   sz: size of reference array e.g. `size(3Darray)`
    %   I: list of linear indices
    %
    % Return:
    %   I1: list of corresponding 1st sub-indices
    %   I2: list of corresponding 2nd sub-indices
    %   I3: list of corresponding 3rd sub-indices
    %

    I1 = rem(I-1, sz(1)) + 1;
    I = (I-I1) / sz(1) + 1;
    I2 = rem(I-1, sz(2)) + 1;
    I3 = (I-I2) / sz(2) + 1;
end