function [source_rays, seg_rays, coords] = ray_cast(interp_source, interp_seg, center, num_rays, ray_interval)
                % * Compute Rays
                % Getting the range limit of the ray which will be the shortest
                % distance from the centre to the bounadry Need to find the
                % limits of the raw
                image_sz = size(interp_source);
                limit_row = abs(image_sz(1) - center(1));
                limit_col = abs(image_sz(2) - center(2));
                ray_length_limit = min(limit_row, limit_col);

                %Compute rays in polar coords
                ray_angle_interval = 2*pi/num_rays;
                radial = 0:ray_interval:ray_length_limit;
                theata = 0:ray_angle_interval:2*pi;

                %Convert rays to cartesian coords
                x_component = radial'*cos(theata) + center(1);
                y_component = radial'*sin(theata) + center(2);

                % * Cast rays
                interp_source = double(interp_source);

                source_rays = interp2(interp_source, x_component(:),...
                    y_component(:));

                seg_rays = interp2(interp_seg, x_component(:),...
                    y_component(:));

                %Need to reshape
                source_rays = ...
                    reshape(source_rays,[size(y_component,1) size(y_component,2)]);
                seg_rays = ...
                    reshape(seg_rays,[size(y_component,1) size(y_component,2)]);
                coords = cat(3, x_component, y_component);
            end