% profiling - run interpolation
if ~exist('AQnet','var')
    run profile_base.m
else
    tic
end

AQnet = AQnet.MakeTubePatches(type='source');
toc