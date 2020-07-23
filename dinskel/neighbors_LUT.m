    function nb = neighbors_LUT(object, conn)
        % 8 connected neighbor positions of every pixel
        % for 2d object, creates a 4d array of size size(object, 1),
        % size(object, 2), 2, 8.
        switch conn
            case 8 % 2D
                nb = zeros([size(object), 2, 8]);
                for i = 1:size(object, 1)
                    for j = 1:size(object, 2)
                        nb(i, j, :, :) = [i-1 j-1; i-1 j; i-1 j+1; ...
                            i j-1; i j+1; i+1 j-1; i+1 j; i+1 j+1;]';
                    end
                end
            case 26 % 3D
                nb = zeros([size(object), 3, 26]);
                for i = 1:size(object, 1)
                    for j = 1:size(object, 2)
                        for k = 1:size(object, 3)
                            % 26 x 3 matrix of all neighbors.
                            nb(i, j, k, :, :) = [i-1 j-1 k-1; ...
                                i-1 j k-1; i-1 j+1 k-1; i j-1 k-1; ...
                                i j k-1; i j+1 k-1; i+1 j-1 k-1; ...
                                i+1 j k-1; i+1 j+1 k-1; i-1 j-1 k; ...
                                i-1 j k; i-1 j+1 k; i j-1 k; ...
                                i j+1 k; i+1 j-1 k; ...
                                i+1 j k; i+1 j+1 k; i-1 j-1 k+1; ...
                                i-1 j k+1; i-1 j+1 k+1; i j-1 k+1; ...
                                i j k+1; i j+1 k+1; i+1 j-1 k+1; ...
                                i+1 j k+1; i+1 j+1 k+1;]';
                        end
                    end
                end
        end
    end