function ijk = RAS2IJK(ras,meta)
    % RAS2IJK converts from RAS to IJK coordinates
    %   ras is a vector (3x1) in RAS coordinates
    %   ijk returned are voxel indices (3x1)
    %   ijk may be a floating point number
    
    % get affine
    rasaff = meta.Transform.T';

    % invert and apply affine
    ijk1 = rasaff \ [ras; 1];
    
    % remove 1
    ijk = ijk1(1:3);

end      