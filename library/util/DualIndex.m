function Out = DualIndex(mat, vec1, vec2)
% By Ashkan Pakzad 2021
% Handles indexing with multiple vectors for a given matrix by utilising
% sub2ind. Can also handle nan index by eliminating in both inputs.

nanvals = isnan(vec1) | isnan(vec2);
vec1 = vec1(~nanvals);
vec2 = vec2(~nanvals);

ind = sub2ind(size(mat), vec1, vec2);

Out = mat(ind);
end
