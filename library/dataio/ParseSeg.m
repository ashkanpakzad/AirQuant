function robustseg = ParseSeg(seg)
    % Fills holes in segmentation and keeps only the largest connected 
    % component.
    %
    % Args:
    %   seg: input segmentation
    % Returns:
    %   robustseg: input segmentation without holes
    
    % fill holes
    segfilled = imfill(seg,'holes');
    
    % preprocess segmentation to keep largest
    % connected component.
    CC = bwconncomp(segfilled);
    numOfPixels = cellfun(@numel,CC.PixelIdxList);
    [~,indexOfMax] = max(numOfPixels);
    robustseg = false(size(segfilled));
    robustseg(CC.PixelIdxList{indexOfMax}) = 1;