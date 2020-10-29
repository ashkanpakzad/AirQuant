function binary_image = ...
    Skeletonise_via_PTK_with_index(binary_image, index_fixed_points)
% This function has been apadted from PTKSkeletonise.

% The input is the binary image and the fixed point- the finxed point will 

% The output is the binary image with a skeletonisation

%% Marking the end point of the image
% Marks endpoints with 3
binary_image = int8(binary_image ~= 0);
binary_image(index_fixed_points) = 3;


%%

binary_image = AddBorder_adapted(binary_image,2);
[i, j, k] = ind2sub([3 3 3], 1 : 27);
direction_vectors = [i' - 2, j' - 2, k' - 2];
raw_image = binary_image;

previous_image = zeros(size(raw_image), 'uint8');

iteration = 0;

while ~isequal(previous_image, raw_image)
    previous_image = raw_image;
    
    iteration = iteration + 1;
    
    
    % For each of the 6 principal directions
    for direction = [5, 23, 11, 17, 13, 15]
        
        direction_vector = direction_vectors(direction,:);
        i = direction_vector(1);
        j = direction_vector(2);
        k = direction_vector(3);
        
        % Detect border points and get their indices
        [b_i, b_j, b_k] = ind2sub(size(raw_image) - [2 2 2], ...
            find((1 == raw_image(2:end-1, 2:end-1, 2:end-1)) & (0 == raw_image(2+i:end-1+i,2+j:end-1+j,2+k:end-1+k))));
        b_i = b_i + 1; b_j = b_j + 1; b_k = b_k + 1;
        
        % Iterate through each border point and delete (set to zero) if
        % it is a simple point
        for ii = 1 : length(b_i)
            i = b_i(ii); j = b_j(ii); k = b_k(ii);
            raw_image(b_i(ii), b_j(ii), b_k(ii)) = ~PTKIsSimplePoint(binary_image(i-1:i+1, j-1:j+1, k-1:k+1));
        end
    end
end

binary_image = ChangeRawImage_adapted(binary_image,raw_image);
% remove border
border_size = 2;
binary_image = binary_image(1+border_size : end-border_size,...
    1+border_size : end-border_size, 1+border_size : end-border_size);

end
