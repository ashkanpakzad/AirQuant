function grid = grid_preview(tar_file, n_rows, n_cols, savepath)
    % Given a tar file of patches, display a grid of them and save as png.
    
    if nargin < 4
        savepath = 0;
    else
        savepath = parse_filename_extension(savepath, '.png');
    end

    % Load the tar file
    atemp_dir = tempname;

    % check tar path
    untar(tar_file,atemp_dir)

    % Get the number of patches
    patches = dir(fullfile(atemp_dir, '*.tif'));
    npatches = length(patches);
    randi = randperm(npatches);

    % check enough patches for grid
    if npatches < n_rows*n_cols
        warning('Not enough patches in tar file for grid. Resizing grid.')
        n_cols = floor(npatches/n_rows);
    end

    % get the size of the patches
    patch_size = size(imread(fullfile(atemp_dir,patches(1).name)));

    % synthesize the grid
    grid = zeros(patch_size(1)*n_rows, patch_size(2)*n_cols);
    k=1;
    for i=1:n_rows
        for j=1:n_cols
            grid((i-1)*patch_size(1)+1:i*patch_size(1), (j-1)*patch_size(2)+1:j*patch_size(2)) = imread(fullfile(atemp_dir,patches(randi(k)).name));
            k=k+1;
        end
    end

    % delete the temp dir
    rmdir(atemp_dir, 's');

    % balance the grid
    grid = imadjust(rescale(grid));

    % save the grid as png
    if savepath
        imwrite(single(grid), savepath);
    end
end