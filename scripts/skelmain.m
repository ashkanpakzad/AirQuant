% By Ashkan Pakzad, 2020. ashkanpakzad.github.io
    % Compute skeleton of input airway graph using PTKskeleton.
    % This method ensures there are no internal loops within the skeleton.
    % Dependent on PTK github
casenames = ["github_demo"]; 
AirQuantAddPath(); % Ensure all AirQuant files are in matlab path
for ii = 1:length(casenames)
    segname = [char(casenames(ii)), '_seg.nii.gz'];
    % must be in matlab current path
    PTKskel(segname, 1);
    
    %% check
    % assumes gzip compression
    skelname = [segname(1:end-7), '_PTKskel.nii.gz'];
    
    seg = niftiread(segname);
    skel = niftiread(skelname);
    
    figure;
    patch(isosurface(seg),'EdgeColor', 'none','FaceAlpha',0.3);
    title(char(casenames(ii)))
    hold on
    isosurface(skel)
    axis vis3d
    view(80,0)
end