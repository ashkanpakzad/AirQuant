function [t_points, s] = spline_points(spline, sample_interval)
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

% Construct arclength interval array s to total arclength
s = 0:sample_interval:L(end);

% compute corresponding parametrised t for s
t_points = zeros(size(s));

for i = 2:length(s)    
    % identify spanned segments
    for j = 2:length(L)
        if s(i-1) >= L(j-1) && s(i-1) < L(j)
            pp_0 = j-1;
        end
        if s(i) >= L(j-1) && s(i) < L(j)
            pp_1 = j-1;
        end
    end
        
    % compute arclength to calculate in final segment
    if pp_0 ~= pp_1
        l = sample_interval - (L(pp_1) - s(i-1));
        t0 = breaks(pp_1);
    else
        l = sample_interval;
        t0 = t_points(i-1);
    end
        
    % setup arclength function
    l_integral = arclengthfunc(pp_1);
    
    % set up integral such that it equals to 0
    l_integral0 = @(t) integral(l_integral,t0,t) - l;
    % identify where t=0 between expected intervals.
    t1 = fzero(l_integral0, [t0, breaks(pp_1+1)]);
    t_points(i) = t1;

end


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