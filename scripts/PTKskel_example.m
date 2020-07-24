% By Ashkan Pakzad, 2020. ashkanpakzad.github.io
    % Compute skeleton of input airway graph using PTKskeleton.
    % This method ensures there are no internal loops within the skeleton.
    % Dependent on PTK github
    
Airwaysegname = 'github_demo_seg.nii.gz';
AirQuantAddPath(); % Ensure all AirQuant files are in matlab path

PTKskel(Airwaysegname);