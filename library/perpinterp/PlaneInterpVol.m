function out_plane = PlaneInterpVol(vol, voxdim, point, normal, options)
    % short desc
    %
    % long desc
    %
    % .. todo: add documentation to this function
    %
    % Args:
    %   vol(`single` or `double`): volume to interpolate.
    %   voxdim(`triplet`): voxel dimensions of :attr:`vol`
    %   point(`triplet`): point at which to define interpolation plane.
    %   normal(`triplet`): normal vector of interpolation plane.
    %   plane_sz: `OPTIONAL, default = 80` pixel size of plane.
    %   sample_sz: `OPTIONAL, default = 0.5` interpolation sample size.
    %   method: `OPTIONAL, default = 'cubic'` interpolation method. see
    %       the same arg in MATLAB's `interp3` for more details.
    %   offgrid_val: `OPTIONAL, default = 0` value to set plane pixels that
    %   are out of volume.
    %
    % Return:
    %   out_plane(type):
    %

    % compute interpolated slice size
    % * Interpolate Perpendicular Slice per spline point
    % Interpolates a CT plane of the image.
    % Based on original function by Kin Quan 2018

    % Construct vol grid
    arguments
        vol 
        voxdim {mustBeNumeric}
        point (3,1) {mustBeNumeric}
        normal (3,1) {mustBeNumeric}
        options.plane_sz (1,1) = 40
        options.sample_sz (1,1) = 0.5
        options.method char = 'cubic'
        options.offgrid_val (1,1) {mustBeNumeric} = 0
    end

    image_sz = single(size(vol));
    [x_domain , y_domain , z_domain] = ...
        meshgrid(1:image_sz(2),1:image_sz(1),1:image_sz(3));
    x_domain = x_domain*voxdim(1);
    y_domain = y_domain*voxdim(2);
    z_domain = z_domain*voxdim(3);

    % Get plane grid
    basis_vecs = Orthonormal_basis_with_tangent_vector(normal);
    plane_grid = Grids_coords_for_plane(basis_vecs(:,3),...
        basis_vecs(:,2), point, options.plane_sz, options.sample_sz);
        
    plane_intensities = AQinterp3(x_domain,y_domain,z_domain,vol, ...
        plane_grid,options.method);
    
    % Reshape to plane
    plane_length = sqrt(length(plane_grid(2,:)));
    out_plane = reshape(plane_intensities,...
        [plane_length plane_length]);
    % Replace NaN entries in image with offgridval.
    out_plane(isnan(out_plane)) = options.offgrid_val;