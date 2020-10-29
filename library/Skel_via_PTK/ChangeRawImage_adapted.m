function output_image = ...
    ChangeRawImage_adapted(current_binary_image , new_image)
%% Adapted function from the PTK code

%The input will be the raw image with the new and old

%The output will be the changed image

%% Need to place extra variables

current_image_size = size(current_binary_image);

%%

% Replaces the underlying raw image data. This function is used for
% applying data which relates to the same original image; hence the
% image size must be the same.

% Compare image sizes, but ignore anything beyond the first 3
% dimensions so that we allow switching between quiver plots
% (vectors) and scalars
new_image_size = size(new_image);

% If the length of the 3rd dimension of 'new_image' is 1, Matlab
% will remove that dimension from the size argument, so we need to
% add it in so we can make a proper comparison
if length(new_image_size) == 2
    new_image_size = [new_image_size 1];
end

if length(new_image_size) > 3
    new_image_size = new_image_size(1:3);
end
this_image_size = current_image_size;
if length(this_image_size) > 3
    this_image_size = this_image_size(1:3);
end

if ~isequal(this_image_size, new_image_size)
    error('The new image data must be the same size as the original data');
end

output_image = new_image;

end