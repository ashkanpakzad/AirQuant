function ptcloud = convert_xyz(image, name)
% convert binary image input into xyz text files for skeletonisation.
% 
% BB.txt
% voxel_size
% x_BB0 y_BB0 z_BB0
% x_BB1 y_BB1 z_BB1
% 
% xyz.txt
% num_elem
% x_1 y_1 z_1
% x_2 y_2 z_2
% ... ... ...
% x_n y_n z_n
% 

% config
voxsize = 1;
BB_name = fullfile(name,'BB.txt');
xyz_name = fullfile(name,'xyz.txt');
ply_name = fullfile(name,name);


% make directory
mkdir(name)

% make BB.txt
BB_line1 = voxsize;
writematrix(BB_line1,BB_name, 'FileType', 'text')

BB_line2 = [0, 0 ,0];
BB_line3 = size(image);
BB_line2_3 = [BB_line2; BB_line3];
writematrix(BB_line2_3,BB_name, 'WriteMode','append', 'Delimiter',' ')


% make 0.txt
zero_line1 = sum(image(:));
writematrix(zero_line1,xyz_name, 'FileType', 'text')

k = find(image);
[X, Y, Z] = ind2sub(size(image),k);
zero_line2_n = [X, Y, Z];
writematrix(zero_line2_n,xyz_name, 'WriteMode','append', 'Delimiter',' ')

% make ply
pt_color = zeros(size(zero_line2_n));
pt_color(:,3) = 255;
ptcloud = pointCloud(zero_line2_n,Color=pt_color);

pcwrite(ptcloud,ply_name,PLYFormat="binary");

end