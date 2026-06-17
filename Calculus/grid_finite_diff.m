function [M, x, r, a] = grid_finite_diff(N, L, order, s, options)
% Creates a finite difference differentiation matrix for derivatives of the
% specified order and a corresponding grid vector.
%
% Inputs:
% - N: number of gridpoints (default: 101)
% - L: grid limits, [L(1) L(2)], if entered as a scalar limits are [0 L] (default: 1)
% - order: order of derivative (default: 1)
% - s: finite difference stencil (default: -order:order);
%           1 - central finite differences, s = [-1 0 1] (O(h^2) for order = 1)
%           2 - forward finite difference, s = [0 1] (O(h) for order = 1)
%           3 - backward finite difference, s = [-1 0] (O(h) for order = 1)
%           s - custom, vector of length S for S > order, e.g. s = [-2 -1 0 1 2]
% - options:
%   - periodic: true - grid is periodic, false - grid is not periodic (default)
%   - sparse: true - output as a sparse matrix , false - output as a full (default)
%   - boundary: method used to calculate deriviatives near boundaries;
%           1 - stencil; automatically calculate a stencil of similar size (default)
%           2 - interpolate; linear interpolation
%           0 - zero; replace boundary rows with 0
%   - flip: 1 - flips the basis and differentiation matrix (i.e. order of points reversed)
%           2 - flips the basis and differentiation matrix if x(end) < x(1)
%           0 - does not flip (default)
%
% Outputs:
% - M: differentiation matrix, size [N N]
% - x: grid, size [N 1]
% - r: indices of replaced rows
% - a: finite difference coefficients
%
% ----------------------------------------------------------------------------
% Note: The order of accuracy is generally O(h^(size(s)-order)). The
%       coefficients may be extracted via output 'a' if theoretical
%       consideration is required.
%
% Note: If the grid is not periodic, the derivatives near the endpoints are
%       determined by either a new stencil (boundary = 1) or linear 
%       interpolation (boundary = 2). The new stencil is chosen to be the
%       same size as s so the result will generally be the same order of
%       accuracy. Linear interpolation may work better for low order
%       accuracy derivatives however will generally restrict the accuracy
%       to O(h).
%
% Note: If using this function to create derivative matrices for solving
%       differential equations, the boundary method may be irrelevant as
%       boundary rows are often replaced by boundary conditions. In this
%       case boundary = 0 may be used to set boundary rows to 0.
%
% Note: Outputting M as sparse (sparse = true) can save a lot of memory and
%       speed up calculations if N is large.
%
% Note: This function requires the additional function 'circulant.m'.
% ----------------------------------------------------------------------------

arguments
    N                (1,1) double {mustBeInteger}                         = 101
    L                (:,1) double                                         = 1
    order            (1,1) double {mustBeInteger}                         = 1
    s                (1,:) double {mustBeInteger,StencilCheck(s,N,order)} = -order:order
    options.periodic (1,1) logical                                        = false
    options.sparse   (1,1) logical                                        = false
    options.boundary (1,1) double {mustBeInteger}                         = 1
    options.flip     (1,1) double {mustBeInteger}                         = 0
end

% Define preset stencils:

if isscalar(s)
    switch s
        case 1
            s = [-1 0 1];
        case 2
            s = [0 1];
        case 3
            s = [-1 0];
        otherwise
            eid = 'Stencil:Undefined';
            msg = 'No preset stencil exists for this value.';
            error(eid, msg)
    end
end

% Convert scalar L to [0 L]:

if length(L) == 1
    L = [0 L];
end

% Create vector x:

if options.periodic
    x = linspace(L(1), L(2) - (L(2) - L(1))/N, N)';
else
    x = linspace(L(1), L(2), N)';
end

% Define grid spacing:

h = x(2) - x(1);

% Create linear system for building differentiation matrix:

exp = (0:length(s)-1)';             % define exponent for calculating matrix system
v = factorial(order)*(order==exp);  % create RHS vector for linear system
a = ((s.^exp) \ v)';                % solve linear system for finite difference coefficients

s = s(abs(a) > eps);                % ignore coefficients of zero
a1 = a(abs(a) > eps);

% Build differentiation matrix:

M = zeros(N);
for i = 1:length(s)
    M = M + a1(i) * circulant(N, s(i));
end

% For non-periodic grids, calculate boundary stencils and rows:

if ~options.periodic

    r = [1:max(0,-s(1)) (N+min(-s(end),0)+1):N];    % indices of rows to replace
    M(r,:) = 0;                                     % set replaced rows to zero

    % Use boundary stencils to calculate derivatives at boundaries:

    if options.boundary == 1

        for ri = r      % loop through replaced rows

            % Create boundary stencils:

            sm = ri + [s(1) s(end)];

            if sm(1) < 1; sb = (s(1):s(end))+1-sm(1); end
            if sm(2) > N; sb = (s(1):s(end))+N-sm(2); end

            % Create finite difference coefficients for boundary stencils:

            exp = (0:length(sb)-1)';
            v = factorial(order)*(order==exp);
            a2 = ((sb .^ exp) \ v)';

            % Set boundary rows of differentiation matrix:

            M(ri, sb + ri) = a2;

        end

    end

    % Use linear interpolation to calculate derivatives at boundaries:

    if options.boundary == 2

        % Interpolate outwards to boundary rows at each domain end:

        for ri = max(0,-s(1)):-1:1
            M(ri,:) = 2 * M(ri+1,:) - M(ri+2,:);
        end

        for ri = (N+min(-s(end),0)+1):N
            M(ri,:) = 2 * M(ri-1,:) - M(ri-2,:);
        end

    end

else

    % No replaced rows for periodic case:

    r = [];

end

% Divide through by factors of h:

M = (1 / h) ^ order * M;

% Reverse rows of x and M if options.flip = 1 or options.flip = 2 and x(end) < x(1):

if options.flip == 1 || (options.flip == 2 && x(end) < x(1))
    x = x(end:-1:1);
    M = M(end:-1:1,end:-1:1);
end

% Convert M to sparse of options.sparse = true:

if options.sparse
    M = sparse(M);
end

end

function StencilCheck(s, N, order)
    if order >= length(s) && ~(isequal(s, 1) && order < 3)
        eid = 'Order:InconsistentStencil';
        msg = 'The number of stencil points must exceed the order of the derivative.';
        error(eid, msg)
    end
    if N < (1 + s(end) - s(1))
        eid = 'Stencil:InconsistentGrid';
        msg = 'The number of stencil points must not exceed the number of grid points.';
        error(eid, msg)
    end

end