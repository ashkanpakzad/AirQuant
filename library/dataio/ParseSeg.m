function seg = ParseSeg(seg, options)
    % Fills holes in segmentation and keeps only the largest connected
    % component.
    %
    % Args:
    %   seg: input segmentation
    % Returns:
    %   robustseg: input segmentation without holes
    
    arguments
    seg
    options.fillholes = 1
    options.largestCC = 1
    end
    
    % fill holes
    if options.fillholes == 1
        seg = imfill(seg,'holes');
    else
        return
    end

    % preprocess segmentation to keep largest
    % connected component.
    if options.largestCC == 1
        CC = bwconncomp(seg);
        numOfPixels = cellfun(@numel,CC.PixelIdxList);
        [~,indexOfMax] = max(numOfPixels);
        seg = false(size(seg));
        seg(CC.PixelIdxList{indexOfMax}) = 1;
    else
        return
    end