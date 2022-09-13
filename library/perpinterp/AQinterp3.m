function plane_intensities = AQinterp3(x_domain,y_domain,z_domain,vol,plane_grid,method,gpu)


if gpu
    assert(strcmp(method,'linear'),'Can only use "linear" method with gpu optimisation.')
    assert(license('test','distrib_computing_toolbox'),' "Parallel Computing Toolbox required for gpu optimisation.')

    % use gpu
    x_domain = gpuArray(x_domain);
    y_domain = gpuArray(y_domain);
    z_domain = gpuArray(z_domain);
    vol = gpuArray(vol);
    plane_grid = gpuArray(plane_grid);

    plane_intensities = interp3(x_domain,y_domain,z_domain,...
        vol,plane_grid(2,:),plane_grid(1,:),...
        plane_grid(3,:),'linear');
    plane_intensities = gather(plane_intensities);
else
    % use standard
    plane_intensities = interp3(x_domain,y_domain,z_domain,...
        vol,plane_grid(2,:),plane_grid(1,:),...
        plane_grid(3,:),method);
end




end