function [min_vol,min_point] = min_tube_vol(vol,point,zooms,plane_size)
    % Returns the minimum volume required t make perpendicular slices of
    % tubes.
    
    arguments
        vol (:,:,:)
        point (3,1)
        zooms (3,1)
        plane_size (1,1)
    end

    % identify indices of bounding box
    min_BB_mm = point - plane_size/2;
    max_BB_mm = point + plane_size/2;
    
    % transform to index
    min_BB = floor(min_BB_mm./zooms);
    max_BB = ceil(max_BB_mm./zooms);
    
    % check BB does not exceed vol extremes
    for ii = 1:length(min_BB)
        if min_BB(ii) < 1
            min_BB(ii) = 1;
        end

        maxdim = size(vol,ii);
        if max_BB(ii) > maxdim
            max_BB(ii) = maxdim;
        end
    end
  
    % transform vol
    min_vol = vol(min_BB(1):max_BB(1),min_BB(2):max_BB(2),...
        min_BB(3):max_BB(3));

    % transform point, plus zooms due to index 1
    min_point = point - min_BB.*zooms+zooms;
end

