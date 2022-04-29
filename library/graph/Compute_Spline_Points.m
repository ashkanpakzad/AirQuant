function varargout = Compute_Spline_Points(spline, sample_interval)
    % Computes the parametrised points along a spline at sample interval.
    % Polynomial splines are defined by discontinuous polynomial
    % functions. They are parametrised by a continuous variable. An
    % interval in parameterised space does not correspond to the same
    % interval in cartesian space. This function identifies the
    % parametrised points for an equidistant spacing along the arclength of
    % the spline.
    %
    % .. note: The last value in :attr:`arcpoints` is equal to 
    % :attr:`total_arclength` if the latter is exactly divisible by  
    % :attr:`sample_interval`. 
    % 
    % Example:
    %   totalarclength = Compute_Spline_Points(spline);
    %   [totalarclength, para_points, arc_points] = Compute_Spline_Points(spline, sample_interval);
    %
    % Args:
    %   spline: polynomial spline object, e.g. computed using `cscvn`.
    %   sample_interval: scalar value defining the sample interval in
    %   cartesian space.
    %
    % Return:
    %   total_arclength: Total arclength of the 
    %   parapoints: parametrised points on spline corresponding to the
    %       defined sample interval.
    %   arcpoints: cartesian points on spline corresponding to the defined
    %       sample interval.
    
    arguments
        spline 
        sample_interval = []
    end

    
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
        L(i) = ComputeArcLength_pp(breaks(i-1), breaks(i)) + L(i-1);
    end
    total_arclength = L(end);
    varargout{1} = total_arclength;


    if ~isempty(sample_interval)

        % Construct arclength interval array s to total arclength
        arc_points = 0:sample_interval:total_arclength;

        % compute corresponding parametrised t for s
        para_points = zeros(size(arc_points));

        for i = 2:length(arc_points)
            % identify spanned segments
            for j = 2:length(L)
                if arc_points(i-1) >= L(j-1) && arc_points(i-1) < L(j)
                    pp_0 = j-1;
                end
                if arc_points(i) >= L(j-1) && arc_points(i) < L(j)
                    pp_1 = j-1;
                end
            end

            % compute arclength to calculate in final segment
            if pp_0 ~= pp_1
                l = sample_interval - (L(pp_1) - arc_points(i-1));
                t0 = breaks(pp_1);
            else
                l = sample_interval;
                t0 = para_points(i-1);
            end

            % setup arclength function
            l_integral = arclengthfunc(pp_1);

            % set up integral such that it equals to 0
            l_integral0 = @(t) integral(l_integral,t0,t) - l;
            % identify where t=0 between expected intervals.
            l_int = t0;
            u_int = breaks(pp_1+1);
            % ensure left and right intervals differ in sign
            while l_integral0(l_int) > 0 || l_integral0(u_int) < 0
                l_int = l_int - 5;
                u_int = u_int + 5;
            end
            t1 = fzero(l_integral0, [l_int, u_int]);
            para_points(i) = t1;

        end
        varargout{2} = para_points;
        varargout{3} = arc_points;
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