function [Int, M] = integ_mat(F, x, dim, F0, method, options)
% Calculates the integral of F dx by integrating F along the specified dimension
%
% Inputs:
% - F: n-dimensional array
% - x: vector of x values, must match length of integrated dimension
% - dim: dimension of F to integrate along
% - F0: initial value of integral I at lower limit, x = x(1) (default: 0)
% - method:
%   - 0: right-hand Riemann sum for non-constant dx
%   - 1: trapezoidal rule for non-constant dx (default)
%   - 2: Fourier transform
%   - 3: higher order finite difference inversion
% - options:
%   - points: number of points for Newton-Cotes formula (default: 3) (method 3)
%   - order: order of integral (default: 1) (method 2 and 3)
%   - periodic: true/false - F is/isn't periodic (default: false) (method 2)
%   - boundary: method used for boundaries, see 'grid_finite_diff.m' (method 3)
%   - definite: true/false - calculate definite integral over x (default: false)
%
% Outputs:
% - Int: integral of F
% - M: order 1 Newton-Cotes integration matrix (method 3 only)
%
% -------------------------------------------------------------------------
% Notes:
%
% Method 1:
% This method reduces to the regular (second order accurate) trapezoidal
% rule for constant dx but will not be second order accurate for
% non-constant dx.
%
% Method 2:
% Method 2 requires constant grid spacing and uses Fourier transforms,
% requiring F to be periodic in the integrated dimension only. It also
% requires that F integrates to zero over the full dimension. If not true,
% this condition will be enforced by setting the zero Fourier mode to zero.
% This method may perform multiple integrals using the 'order' parameter.
% Integrals of all orders will also have the property that have a domain
% integral of zero. A non-zero value of F0 may be added to give a final
% result for I which does not integrate to zero. However, only one
% integration constant is applied as applying integration constants at each
% integration would lead to linear, quadratic etc terms which are
% non-periodic. Using method = 2 requires 'FFT_grid.m', 'FFT_forward.m' and
% 'FFT_inverse.m'.
%
% Method 3:
% Method 3 require constant grid spacing and uses a Newton-Cotes formula to
% calculate integrals. The number of points used in the Newton-Cotes
% formula is set by the 'points' parameter. Setting points = 2 is
% equivalent to method 1. By default, method 3 uses a Simpson's rules
% result (points = 3). Integral values at odd points may therefore be more
% accurate than those at even points which also use a single trapezoidal
% rule step. Method 3 may be used to do multiple integrals using the 'order'
% parameter. The conditions on I are set using F0 where F0(1, :) is I(1, :)
% and F0(i+1, :) is d^i/dx^i I(1, :) for i = 1:order-1. Using method 3
% requires 'NewtonCotes_coeff.m'.
%
% Note: The definite integral, i.e the integral of F between x(1) and
% x(end), may be calculated by setting definite = true.
% -------------------------------------------------------------------------

arguments
    F                      double
    x                (:,1) double                                        = 1:size(F,1)
    dim              (1,1) double {mustBeInteger,ConsistentDim(F,x,dim)} = 1
    F0                     double                                        = 0
    method           (1,1) double {mustBeInteger,IsConstdx(x,method)}    = 1
    options.points   (:,1) double                                        = 3
    options.order    (1,1) double {mustBeInteger}                        = 1
    options.periodic (1,1) logical                                       = false
    options.boundary (1,1) double {mustBeInteger}                        = 1
    options.definite (1,1) logical                                       = false
end

% Permute F to have differentiated dimension in the first dimension:

N = length(size(F));        % number of dimensions

I = 1:N; I(dim) = [];       % index of non-differentiated dimensions
F = permute(F, [dim I]);
s = size(F);

% Reshape F to a 2D array with all non-differentiated dims in dim 2:

F = reshape(F, s(1), []);

% Permute, check and reshape integration constant array:

if numel(F0) == 1
    F0 = F0 * ones(1, numel(F(1,:)));
end

F0 = permute(F0, [dim I]);
s1 = size(F0, 1);

if (s1 == 1 && method ~= 3) || (s1 == options.order && method == 3)
    F0 = reshape(F0, s1, []);
else
    eid = 'IntegConstant:Inconsistent';
    msg = ['The size of F0 is inconsistent with the ' ...
        'size of F or with the order of the integral (for method = 3).'];
    error(eid, msg)
end

% Integrate along the first dimension:

switch method

    case 0 % simple (right-hand) Riemann sum

        Int = zeros(size(F));

        for i = 2:size(F, 1)
            Int(i, :) = Int(i - 1, :) + (x(i) - x(i-1)) * F(i, :);
        end

        Int = Int + F0;

    case 1 % trapezoidal rule

        Int = zeros(size(F));

        for i = 2:size(F, 1)
            Int(i, :) = Int(i - 1, :) + (x(i) - x(i-1)) * (F(i-1, :) + F(i, :)) / 2;
        end

        Int = Int + F0;

    case 2 % fast Fourier transform

        if ~options.periodic
            eid = 'Grid:NonPeriodic';
            msg = 'Direction must be periodic for this method.';
            error(eid, msg)
        end

        % Define number of grid points and grid limits:

        N = length(x);
        L = [x(1) x(end)+mean(diff(x))];

        [~, k] = FFT_grid(N, L);

        % Use FFTs to calculate integral:

        kn_inv = (1i * k').^-options.order;
        kn_inv(isinf(kn_inv)) = 0;
        Int = real(FFT_inverse(kn_inv .* FFT_forward(F, 1), 1)) + F0;

    case 3 % n-point extended Newton-Cotes formula

        % Build Newton-Cotes integration matrix

        N = length(x);
        h = mean(diff(x));

        M = zeros(N);

        for i = 2:options.points-1
            M(i, 1:i) = NewtonCotes_coeff(i);
        end

        H = NewtonCotes_coeff(options.points);

        for i = options.points:N
            M(i, :) = [zeros(1, i - options.points) H zeros(1, N - i)] + ...
                M(i - options.points + 1, :);
        end

        % Iteratively integrate by applying M then adding the integration constant:

        Int = F;

        for i = 1:options.order
            Int = h * M * Int + F0(options.order - i + 1, :);
        end

    otherwise   % unimplemented methods

        eid = 'Method:NotImplemented';
        msg = ['Using method = ' num2str(method) ' has not been implemented.'];
        error(eid, msg)

end

% If calculating a definite integral, take only the value of I at x(end):

if options.definite
    if options.periodic
        Int = Int(1, :);
    else
        Int = Int(end, :);
    end
end

% Reshape and permute dF to match original F:

Int = reshape(Int, [size(Int, 1) s(2:end)]);
Int = permute(Int, [I(I<dim)+1 1 I(I>dim)]);

% Define M if method ~= 3:

if method ~= 3
    M = [];
end

end

function ConsistentDim(F, x, dim)
    s = size(F);
    if length(x) ~= s(dim)
        eid = 'Size:Mismatch';
        msg = ['Length of x must match length of dimension ' num2str(dim) ' of F'];
        error(eid, msg)
    end
end

function IsConstdx(x, method)
    if max(abs(diff(diff(x)))) > 1e-10 * max(abs(diff(x))) && ismember(method, [2 3])
        eid = 'Spacing:NonConstant';
        msg = ['Vector x must have constant spacing for method = ' num2str(method)];
        error(eid, msg)
    end
end