function [a, b, f_pade, err] = Pade(f, x, options)
% Calculates the Pade approximant of order (n, m) to f(x).
%
% Inputs:
% - f: function entered as function handle or column vector
% - x: vector describing x interval or column vector corresponding to f (default: 1:length(f))
% - options:
%   - m: order of numerator (default: 2)
%   - n: order of denominator (default: 2)
%   - num_exp: row vector of exponents in numerator (default: 0:m)
%   - den_exp: row vector of exponents of non-constant terms in denominator (default: 1:n)
%   - debug: if true, debug information is shown
%   - round: number of significant figures for rounding coefficients, if 0 (default), no rounding 
%   - calc_x: if true, x is calculated using x(1) and x(end) as the endpoints (default: false)
%   - num_pts: number of x points to use if calc_x = true, default set so problem is square
%   - spacing: spacing type if using 'calc_x = true', 'lin' (default) or 'log'
%   - interp: interpolation method (see 'interp1') if calc_x = true (default: 'spline')
%
% Outputs:
% - a: numerator coefficients
% - b: denominator coefficients (including leading constant coefficient 1)
% - f_pade: Pade approximant as a function handle
% - err: error in Pade approximant evaluated over x (sqrt(MSE))
%
% -------------------------------------------------------------------------
% Notes:
% num_exp and den_exp can be explicitly set to control what powers of x
% appear in the numerator and denominator of the Pade approximant. There is
% always a constant term (1) in the denominator. If these are not set, the
% powers of x will be determined from the order of the numerator (m) and
% denominator (n).
%
% Setting calc_x = true can be used to calculate x rather than using the x 
% vector entered. If calc_x = true, x will be determined over the interval 
% [x(1) x(end)] with a number of points set by 'num_pts'. If 'num_pts' is
% not specified or 0, it will be set such that the problem is square, i.e.
% num_pts = length(num_exp) + length(den_exp). The x points are distributed
% linearly if spacing = 'lin', and logarithmically if spacing 'log'. If the
% problem is square, the coefficients can be solved for exactly, while if
% non-squre they will be determined by a least-squares approach. If calc_x
% = true, then f will be determined by direct evaluation if f is a function
% handle, or by interpolation if f is a vector. The interpolation method is
% set using 'interp', with accepted methods given in the MATLAB 'interp1'
% documentation.
%
% Setting round = 0 performs no rounding of the coefficients, a value of
% round > 0, rounds the coefficients to that number of significant figures.
%
% This function requires 'sig_fig.m'.
% -------------------------------------------------------------------------

arguments
    f               (:,1)
    x               (:,1) double {ConsistentType(x,f)} = 1:length(f)
    options.m       (1,1) double {mustBeInteger}       = 2
    options.n       (1,1) double {mustBeInteger}       = 2
    options.num_exp (1,:) double                       = NaN
    options.den_exp (1,:) double                       = NaN
    options.debug   (1,1) logical                      = false
    options.round   (1,1) double {mustBeInteger}       = 0
    options.calc_x  (1,1) logical                      = false
    options.num_pts (1,1) double {mustBeInteger}       = 0
    options.spacing (1,:) char                         = 'lin' % or 'log'
    options.interp  (1,:) char                         = 'spline'
end

% set numerator exponents if not specified:
if sum(isnan(options.num_exp))
    num_exp = 0:options.m;
else
    num_exp = options.num_exp;
end

% set denominator exponents if not specified:
if sum(isnan(options.den_exp))
    den_exp = 1:options.n;
else
    den_exp = options.den_exp;
end

% if calculating x from a given interval, determine x and f:
if options.calc_x

    % set number of points to linear problem size if not specified:
    if options.num_pts == 0
        num_pts = length(num_exp) + length(den_exp);
    else
        num_pts = options.num_pts;
    end

    % distribute points linearly or logarithmically:
    switch options.spacing
        case 'lin'
            xp = linspace(x(1), x(end), num_pts)';
        case 'log'
            if min(x) <= 0
                eid = 'Spacing:NegativeValues';
                msg = 'Non-positive values of x are not permitted for logarithmic spacing.';
                error(eid, msg)
            end
            xp = logspace(log(x(1))/log(10), log(x(end))/log(10), num_pts)';
        otherwise
            eid = 'Spacing:UnrecognisedType';
            msg = 'Value of ''spacing'' not recognised.';
            error(eid, msg)
    end

    % evaluate f values on x using direct evaluation or interpolation:
    if isa(f, "function_handle")
        f = f(xp);
    else
        f = interp1(x, f, xp, options.interp);
    end

    x = xp;

end

% convert function inputs to vectors:
if isa(f, "function_handle")
    f = f(x);
end

% debug information on problem:
if options.debug
    disp([' - grid length: ' num2str(length(x))])
    disp([' - numerator exponents: ' num2str(num_exp)])
    disp([' - denominator exponents: ' num2str([0 den_exp])])
end

% calculate linear problem and solve for coefficients:
M = [x.^num_exp -f.*x.^den_exp];
C = M \ f;

% round coefficients:
if options.round > 0
    for i = 1:length(C)
        C(i) = sig_fig(C(i), options.round);
    end
end

% extract numerator (a) and denominator (b) coefficients:
a = C(1:length(num_exp));
b = [1; C((length(num_exp)+1):end)];

if options.debug
    disp([' - linear problem size: ' num2str(size(M, 1)) ' x ' num2str(size(M, 2))])
end

% determine pade aproximant function:
if nargout > 2
    f_pade = @(x) ((reshape(x, [length(x) 1]) .^ num_exp) * a) ./ ...
        ((reshape(x, [length(x) 1]) .^ [0 den_exp]) * b);

    % calculate error in Pade approximate:
    if nargout > 3
        err = sqrt(sum((f - f_pade(x)).^2 / length(x)));
    end

    % use symbolic expression to simplify f_pade:
    syms x
    f_pade = matlabFunction(simplify(f_pade(x)));
end

end

% validation function for size of f and x:
function ConsistentType(x, f)
    if isa(f, "double") && length(x) ~= length(f)
        eid = 'Length:NotEqual';
        msg = 'Vectors x and f must have the same length.';
        error(eid, msg)
    end
    if isa(f, "function_handle") && length(x) < 2
        eid = 'Length:InsufficientPoints';
        msg = 'Vector x must contain at least two points.';
        error(eid, msg)
    end
end