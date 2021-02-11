function new_image = AddBorder_adapted(binary_image, border_size)
%This is a internal function that has been adepted from object orentated to
%a spcrited based system

%% Need to introduce new variables
image_size_vector = size(binary_image);

%%
%Adds a blank border of border_size voxels to the image in all dimensions
if numel(border_size) == 3
    added_size = border_size;
else
    added_size = [border_size border_size border_size];
end

if islogical(binary_image)
    new_image = false(image_size_vector + 2*added_size);
else
    new_image = zeros(image_size_vector + 2*added_size);
end

new_image(1+added_size(1):end-added_size(1), 1+added_size(2):end-added_size(2), 1+added_size(3):end-added_size(3)) = binary_image;
binary_image = new_image;
%binary_image.CheckForZeroImageSize;
end
%binary_image.Origin = binary_image.Origin - added_size;
%binary_image.NotifyImageChanged;
