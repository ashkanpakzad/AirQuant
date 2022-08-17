function per = ellipsePerimeter(a,b,varargin)
% ELLIPSEPERIMETER Compute the perimeter of an ellipse.
%   per = ELLIPSEPERIMETER(a,b) Computes the perimeter PER of an ellipse 
%   with semi-axes A and B with an infinite series solution set up until 
%   order 20. 
%   per = ELLIPSEPERIMETER(a,b,order) Computes the perimter PER with
%   semi-axes A and B with an infinite series solution set up until the
%   ORDER selected.
%
%   ######################################################################
%   EXAMPLES:
%   >> per = ellipsePerimeter(1,1);
%   Gives as a result: per = 6.283185307179586 = 2*pi
%
%   >> per = ellipsePerimeter(1,5)
%   per = 21.010044539677985 
%    
%   >> per = ellipsePerimeter(1,5,1) 
%   per = 20.943951023931955
%
%   >> per = ellipsePerimeter(1,5,100) 
%   per = 21.010044539689000
%
%   >> r1 = [1 1 1];
%   >> r2 = [4 5 6];
%   >> per = ellipsePerimeter(r1,r2);
%   per = 17.156843550313564   21.010044539677985   24.900079589772904
%
%   >> r2 = [4 5 6];
%   >> per = ellipsePerimeter(1,r2);
%   per = 17.156843550313564   21.010044539677985   24.900079589772904
%
%   ######################################################################
%   The inifite series is given by: 
%               PER = pi(a+b)*SUM(nchoosek(0.5,i)^2 * h^i)              (1)
%   Since nchoosek(n,k) in MATLAB is not defined for non-integers, the
%   definition of the binomial coefficient (and therefore the definition 
%   of the factorial) was changed to use the Gamma function instead.
%                       factorial(n) = gamma(n+1)                       (2)
%   Finally, the binomial coefficient:
%             nchoosek(n,k) = (gamma(n+1)/gamma(n-k+1)/gamma(k+1))      (3)
%   
%   contact: santiago.benito@rub.de

% Parse inputs
p = inputParser;

validationFcn = @(x) isnumeric(x) && sum(x > 0)/size(x,2) == 1;
addRequired(p,'a',validationFcn);
addRequired(p,'b',validationFcn);

defaultOrderValue = 20;
addOptional(p,'order',defaultOrderValue,validationFcn);

parse(p,a,b,varargin{:});

order = p.Results.order;

% Program start
h = (a-b).^2./(a+b).^2;
s = 1;
for i=1:order
    s = s + (h.^i * (gamma(0.5+1)/gamma(0.5-i+1)/gamma(i+1))^2);
end
per = pi*(a+b).*s;
