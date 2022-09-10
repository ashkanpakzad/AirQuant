function angle = vector_angle(P1,P2,units)
% Compute angle of two vectors in either degrees or radians.
%
% Args:
%   P1(vector) = a vector
%   P2(vector) = a vector
%   units(char) = OPTIONAL default = `radians`. Either `radians` or
%     `degrees`.
%
% Returns:
%   1 variable.
%   * angle(`scalar`) = angle between the two input degrees in chosen
%     units.
%
if nargin < 3
    units = 'radians';
end

% compute norms and dots
Y = norm(cross(P1,P2));
X = dot(P1,P2);

switch units
    case 'radians'
    % Angle in radians
    angle = atan2(Y,X); 
    case 'degrees'
    % Angle in degrees
    angle = atan2d(Y,X); 
end