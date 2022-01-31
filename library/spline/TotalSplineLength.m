function totallength = TotalSplineLength(spline)
% Get t break points
breaks = ppbrk(spline, 'breaks');
Npp = ppbrk(spline, 'pieces');
dim = ppbrk(spline, 'dim');

% differentiate spline
diff_spline = fnder(spline, 1);
dcoeff = ppbrk(diff_spline, 'coeff');

% Compute arclength for each segment
L = zeros(size(breaks));
for i = 2:Npp+1
    L(i) = ComputeArcLength_pp(breaks(i-1), breaks(i)) ...
        + L(i-1);
end

totallength = L(end);

function arclength = ComputeArcLength_pp(t0, t1)
% compute the arclength between t0 and t1 where t0 and t1
% are on the same spline segment.

% get spline segment.
% TODO: throw error that t0 and t1 not on same segment if
% fail at this point.
segment = [];
for ii = 1:length(breaks)
    if t0 >= breaks(ii) && t0 < breaks(ii+1)
        segment = ii;
    end
end
% check that t0 and t1 are on the same spline segment
assert(t1 > breaks(segment) && t1<= breaks(segment+1),...
    't0 and t1 need to be on the same spline segment')

I = arclengthfunc(segment);

% compute integral over range
arclength = integral(I,t0,t1);
end

function arclength_integral = arclengthfunc(segment_index)

% Get pp coeff
coeff_index = [1:dim] + dim * (segment_index - 1);
coeff_pp = dcoeff(coeff_index,:);

% construct arc length equation
dd = @(t) 0;
for ii = 1:dim
    current_coeff = coeff_pp(ii,:);
    dd = @(t)(polyval(current_coeff, t - breaks(segment_index))).^2 + dd(t);
end
arclength_integral = @(t) sqrt(dd(t));
end
end
