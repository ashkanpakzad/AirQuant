function [norm_vec, point_3D] = spline_normal(spline, para_point)
    % short desc
    %
    % long desc
    %
    % .. todo: add documentation to this function
    %
    % Args:
    %   x(type):
    %
    % Return:
    %   y(type):
    %
    % Based on original function by Kin Quan 2018
    % interperate real point on spline
    point_3D = fnval(spline, para_point);
    % get tangent of point along spline
    % differentiate along spline to get gradient
    spline_1diff = fnder(spline,1);
    tangent_vec = fnval(spline_1diff, para_point);
    norm_vec = tangent_vec/norm(tangent_vec, 2);
end