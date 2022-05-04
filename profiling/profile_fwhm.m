% profiling - run FWHM measure

if ~exist('AQnet','var')
    run profile_interpsource.m
    run profile_interpseg.m
else
    tic
end

num_rays = 180;
ray_interval = 0.2;
AQnet.Measure('AirwayFWHMesl', num_rays, ray_interval);
toc