function [dF, M] = diff_mat(F, x, dim, method, options)
% Calculates dF/dx by differentiating F along the specified dimension
%
% Inputs:
% - F: n-dimensional array
% - x: vector of x values, must match length of differentiated dimension
% - dim: dimension of F to differentiate along (default: 1)
% - method:
%   - 0: approximate differences (diff(F)/diff(x))
%   - 1: centred finite difference for non-constant dx (default)
%   - 2: Fourier transform
%   - 3: higher order finite difference
% - options:
%   - stencil: finite difference stencil, see 'grid_finite_diff.m' (method 3)
%   - order: order of derivative (default: 1) (method 2 and 3)
%   - periodic: true/false - F is/isn't periodic (default: false) (method 2)
%   - boundary: method used for boundaries, see 'grid_finite_diff.m' (method 3)
%
% Outputs:
% - dF: derivative of F
% - M: order 'order' differentiation matrix (method 3 only)
%
% -------------------------------------------------------------------------
% Notes: methods 2 and 3 require constant grid spacing. Method 1 reduces to
% regular (second order) centred finite differences for constant dx but
% will not be second order accurate for non-constant dx. Method 2 using
% Fourier transforms so requires F to be periodic in the differentiated
% dimension only. Methods 2 and 3 uses additional options and may be used
% to calculate higher order derivatives. Using method = 2 requires
% 'FFT_grid.m', 'FFT_forward.m' and 'FFT_inverse.m' and using method = 3
% requires 'grid_finite_diff.m' and 'circulant.m'.
% -------------------------------------------------------------------------

arguments
    F                      double
    x                (:,1) double                                        = 1:size(F,1)
    dim              (1,1) double {mustBeInteger,ConsistentDim(F,x,dim)} = 1
    method           (1,1) double {mustBeInteger,IsConstdx(x,method)}    = 1
    options.stencil  (:,1) double                                        = [-1 0 1]
    options.order    (1,1) double {mustBeInteger}                        = 1
    options.periodic (1,1) logical                                       = false
    options.boundary (1,1) double {mustBeInteger}                        = 1
end

% Permute F to have differentiated dimension in the first dimension:

N = length(size(F));        % number of dimensions

I = 1:N; I(dim) = [];       % index of non-differentiated dimensions
F = permute(F, [dim I]);
s = size(F);

% Reshape F to a 2D array with all non-differentiated dims in dim 2:

F = reshape(F, s(1), []);

% Differentiate along the first dimension:

switch method

    case 0  % simple approximate finite difference using Matlab 'diff'

        dF = diff(F) ./ diff(x);

    case 1  % centred finite differences with interpolation

        dF = zeros(size(F));

        % Define grid and field differences:

        dx1 = x(2:end-1) - x(1:end-2);
        dx2 = x(3:end)   - x(2:end-1);
        dx3 = x(3:end)   - x(1:end-2);

        dF2 = F(3:end, :)   - F(2:end-1, :);
        dF1 = F(2:end-1, :) - F(1:end-2, :);

        % Calculate interpolated derivative in the domain interior:

        dF(2:end-1, :) = ((dx1 ./ dx2) .* dF2 + (dx2 ./ dx1) .* dF1) ./ dx3;

        % Calculate interpolated derivative at boundaries:

        dF(1, :) = (-dx1(1) * dF2(1, :) / dx2(1) + (dx1(1) + dx3(1)) * dF1(1, :) / dx1(1)) / dx3(1);
        dF(end, :) = ((dx2(end) + dx3(end)) * dF2(end, :) / dx2(end) - dx2(end) * dF1(end, :) / dx1(end)) / dx3(end);

    case 2  % fast Fourier transform

        if ~options.periodic
            eid = 'Grid:NonPeriodic';
            msg = 'Direction must be periodic for this method.';
            error(eid, msg)
        end

        % Define number of grid points and grid limits:

        N = length(x);
        L = [x(1) x(end)+mean(diff(x))];

        [~, k] = FFT_grid(N, L);

        % Use FFTs to calculate derivative:

        dF = real(FFT_inverse((1i * k').^options.order .* FFT_forward(F, 1), 1));

    case 3  % higer order finite difference stencil

        % Define number of grid points and grid limits:

        N = length(x);

        if options.periodic
            L = [x(1) x(end)+mean(diff(x))];
        else
            L = [x(1) x(end)];
        end

        % Create differentiation matrix:

        M = grid_finite_diff(N, L, options.order, options.stencil, ...
            boundary = options.boundary, ...
            periodic = options.periodic, ...
            sparse = true);

        % Apply differentiation matrix to F:

        dF = M * F;

    otherwise   % unimplemented methods

        eid = 'Method:NotImplemented';
        msg = ['Using method = ' num2str(method) ' has not been implemented.'];
        error(eid, msg)

end

% Reshape and permute dF to match original F:

dF = reshape(dF, [size(dF, 1) s(2:end)]);
dF = permute(dF, [I(I<dim)+1 1 I(I>dim)]);

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