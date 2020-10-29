function binary_image = ...
    Skeletonise_via_PTK_with_index(binary_image, index_fixed_points)
% This function has been apadted from PTKSkeletonise.

% The input is the binary image and the fixed point- the finxed point will 

% The output is the binary image with a skeletonisation

%% Marking the end point of the image
% Marks endpoints with 3
binary_image = MarkEndpoints(binary_image, index_fixed_points);


%%

binary_image = AddBorder_adapted(binary_image,2);
direction_vectors = CalculateDirectionVectors;

raw_image = binary_image;

previous_image = zeros(size(raw_image), 'uint8');

iteration = 0;

%Need to set the flag
use_mex_simple_point = false;

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
        for i = 1 : length(b_i)
            raw_image(b_i(i), b_j(i), b_k(i)) = ~IsPointSimple(raw_image, b_i(i), b_j(i), b_k(i), use_mex_simple_point);
            
        end
    end
end
binary_image = ChangeRawImage_adapted(binary_image,raw_image);
binary_image = RemoveBorder_adapted(binary_image,2);
end
