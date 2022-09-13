% profiling - run FWHM measure

tic
num_rays = 180;
ray_interval = 0.2;
AQnet.Measure('AirwayFWHMesl', num_rays, ray_interval);
toc