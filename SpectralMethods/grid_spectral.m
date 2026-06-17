function [M, x] = grid_spectral(options)
% Creates a spectral differentiation matrix, M, and the corresponding 
% vector of grid-points, x. Can also construct grid consisting of multiple
% segments with different point distributions in each.
%
% - N: number of gridpoints in each segment, scalar or vector, default: 16
% - L: interval length or segment endpoints, scalar or vector, default: [0 1]
% - type: character array or cell, accepted types: 'fourier' (default), 
%       'chebyshev1', 'chebyshev2', 'legendre', 'laguerre', 'hermite',
%       'sinc', 'periodicsinc', 'custom'
% - decay: decay scale, Laguerre or Hermite sections only, replaces second
%       segment endpoint if set to finite number, default: NaN
% - flip: if true, segment is flipped before x points are calculated then
%       flipped back, use to reverse asymmetric points such as Laguerre
% - x, w, wp: gridpoints, weight function and weight function derivative
%       for setting custom collocation grids
% - order: order of derivative for Fourier points, default: 1
% - boundary_weight: weighting applied to left and right hand derivatives
%       when setting derivative as segment boundaries, default: [0.5 0.5]
% - debug: if true, function displays debug information about grid
%       construction and parameters
% -------------------------------------------------------------------------
% Grid Types:
%
% Fourier grid:
%
% This method calculates derivative matrix using DFT matrices. Higher order
% derivatives can be specified by setting the 'order' argument. Functions
% defined on this grid should be periodic.
%
% Chebyshev grid of the 1st and 2nd kind:
%
% These methods define Chebyshev points using the cosine function then
% calculate the derivative matrix using the weight function w(x) = 1.
% Chebyshev points of the 2nd kind include the domain endpoints.
%
% Legendre grid:
%
% This method defines Legendre points by determining the eigenvalues of a
% symmetric Jacobi matrix and calculates derivatives using the weight
% function w(x) = 1. These points do not include the endpoints.
%
% Laguerre grid:
%
% This method defines Laguerre points using the eigenvalues of a Jacobi
% matrix and calculates derivatives using the weight function w(x) =
% exp(-x/D) for a decay scale D. D may be set explicitly using the decay
% parameter or will be set automatically using L if D is not specified.
% These points are asymmetric so flip can be used to cluster points on the
% right-hand boundary.
%
% Hermite grid:
%
% This method defines Hermite points using the eigenvalues of a Jacobi
% matrix and calculates derivatives using the weight function w(x) =
% exp(-x^2/D^2) for decay scale D. D may be set explicitly using the decay
% parameter or will be set automatically using L if D is not specified.
%
% Sinc:
%
% This method defines the grid using the semi-discrete Fourier transform by
% representing a function as a weighted sum of sinc functions. Unlike many
% spectral grids, this grid has evenly spaced points. Even N is strongly
% recommended.
%
% Periodic Sinc:
%
% This method uses the same periodic grid as the Fourier method and
% represents functions as a weighted sum of periodic sinc functions. The
% grid and differentiation are identical to the Fourier case. This method
% is only included for demonstration and completeness.
%
% Custom:
%
% This method allows you to calculate a custom spectral collocation matrix
% by specifying the grid, x, weight function evaluated on the grid, w(x),
% and weight function derivative evaluated on the grid, w'(x). Derivatives
% will not be accurate unless the points are clustered near boundaries with
% average spacing O(N^-2) and in the interior with average spacing O(N^-1).
% -------------------------------------------------------------------------
% Notes:
%
% - The method relies on the expansion f(x) ~ Sum_n [f(x_n) phi_n(x) w(x)]
%   therefore the weight function w(x) is the square root of the weight
%   function in the orthogonality relation for the given polynomials.
%
% - Higher order derivative matrices can be calculated similarly however
%   using d^n/dx^n = M^n is generally correct to a high degree of accuracy.
%
% - This script can build composite grids by joining grid segments together
%   provided the endpoints of each segment coincide. Use vector entries for
%   N and L and multi-element cell arrays for type. E.g.
%   grid_spectral(type = {'chebyshev2', 'laguerre'}, N = [4 5], L = [0 1 2])
% -------------------------------------------------------------------------

arguments
    options.N               (1,:) double {isinteger}  = 16
    options.L               (1,:) double              = [0 1]
    options.type            (1,:)        {cellorchar} = 'fourier'
    options.decay           (1,:) double              = NaN
    options.flip            (1,:) logical             = false
    options.x               (:,1) double              = linspace(0, 1, 16)'
    options.w               (:,1) double              = ones(16, 1)
    options.wp              (:,1) double              = zeros(16, 1)
    options.order           (1,:) double              = 1
    options.boundary_weight (1,2) double              = [0.5 0.5]
    options.debug           (1,1) logical             = false
end

% Set number of segments:
S = numel(options.N);

% Extend L to full interval or identify inconsistent length:
if numel(options.L) ~= S + 1
    if numel(options.L) == S && options.L(1) ~= 0
        options.L = [0 options.L];
    else
        error('Size:InconsistentEntry', ['The argument ''L'' should have length ' ...
            num2str(S + 1) ' for length(N) = ' num2str(S)])
    end
end

% Convert type to a cellstr if not already and check length:
type = cellstr(options.type);
if length(type) ~= S
    error('Size:InconsistentEntry', ['The argument ''type'' should have length ' ...
        num2str(S) ' for length(N) = ' num2str(S)])
end

% Set decay vector if scalar specified (same decay for all sections):
if length(options.decay) < S
    decay = options.decay(1) * ones(1, S);
else
    decay = options.decay;
end

% Set decay vector if scalar specified (same decay for all segments):
if length(options.flip) < S
    flip = options.flip(1) * ones(1, S);
else
    flip = options.flip;
end

% Set decay vector if scalar specified (same decay for all segments):
if length(options.order) < S
    order = options.order(1) * ones(1, S);
else
    order = options.order;
end

% Define left and right derivative weightings at interior boundaries:
b1 = options.boundary_weight(1);
b2 = options.boundary_weight(2);

if options.debug
    if S == 1
        disp('Single Segment Grid:')
    else
        disp('Multi-Segment Grid:')
        disp([' - Boundary derivative weights: [' num2str(b1) ', ' num2str(b2) ']'])
    end
end

% Loop through all grid segments:
for n = 1:S

    if options.debug
        if S > 1; disp(['Segment ' num2str(n) ':']); end
        disp([' - Grid points: ' num2str(options.N(n))])
        disp([' - Grid limits: [' num2str(options.L(n)) ', ' num2str(options.L(n+1)) ']'])
    end

    [M_n, x_n] = create_Mx(options.N(n), options.L(n:n+1), type{n}, decay(n), ...
        flip(n), options.x, options.w, options.wp, order(n), options.debug);

    if n == 1

        x = x_n;
        M = M_n;

    else

        if abs(x(end) - x_n(1)) < 10 * eps

            M = [M(1:end-1, :) zeros(length(x)-1, length(x_n)-1); ...
                b1*M(end, 1:end-1) b1*M(end, end)+b2*M_n(1, 1) b2*M_n(1, 2:end); ...
                zeros(length(x_n)-1, length(x)-1) M_n(2:end, :)];

            x = [x; x_n(2:end)];

        else

            error('Endpoints:Mismatch', ['The sub-grid endpoints do not match at x = ' num2str(x(end))])

        end
    end

end

end

% Validation function for type input:
function cellorchar(txt)
    if ~isa(txt, "char") && ~isa(txt, "cell")
        eid = 'Type:Inconsistent';
        msg = 'Argument ''type'' must be a character array or cell.';
        error(eid, msg)
    end
end

% Functin to build M and x for each segment:
function [M, x] = create_Mx(N, L, type, decay, flip, x, w, wp, order, debug)

arguments
    N     (1,1) double
    L     (2,1) double
    type  (1,:) char
    decay (1,1) double
    flip  (1,1) logical
    x     (:,1) double
    w     (:,1) double
    wp    (:,1) double
    order (:,1) double
    debug (1,1) logical
end

% flip L if flip is true:
if flip
    L = L([2 1]);
end

% Create grid by type:
switch type
    
    case 'fourier'
    
        if debug; disp(' - Grid type: Fourier'); end

        Lx = L(2) - L(1);
        dx = Lx / N;
        x = (L(1):dx:L(2)-dx)';
        F = exp(-2 * pi * 1i * (0:N-1)' * (0:N-1) / N); % Fourier transform matrix
        k = 2 * pi / Lx * [0:ceil(N/2-1) -floor(N/2):-1];
        M = real(1 / N * F' * diag(1i*k) .^ order * F);

    case 'sinc'

        if debug; disp(' - Grid type: Sinc'); end

        x = linspace(L(1), L(2), N)';
        dx = x(2) - x(1);
        M = real((-1) .^ ((x - x') / dx) ./ (x - x'));
        M(eye(N) == 1) = 0;

    case 'periodicsinc'

        if debug; disp(' - Grid type: Periodic Sinc'); end

        if rem(N, 2) == 1
            fprintf(2, 'Warning: even values of n are recommended\n')
        end

        Lx = L(2) - L(1);
        dx = Lx / N;
        x = (L(1):dx:L(2)-dx)';
        M = real(pi / Lx * (-1) .^ ((x - x') / dx) ./ tan(pi * (x - x') / Lx));
        M(eye(N) == 1) = 0;

    case {'chebyshev1', 'chebyshev2', 'legendre', 'laguerre', 'hermite', 'custom'}

        switch type
            case 'chebyshev1'

                if debug; disp(' - Grid type: Chebyshev 1'); end

                x = (L(1) + L(2)) / 2 - (L(2) - L(1)) / 2 * cos((1/2:(N-1/2))*pi/N)';
                w = ones(N, 1);
                wp = zeros(N, 1);

            case 'chebyshev2'

                if debug; disp(' - Grid type: Chebyshev 2'); end

                x = (L(1) + L(2)) / 2 - (L(2) - L(1)) / 2 * cos((0:(N-1))*pi/(N-1))';
                w = ones(N, 1);
                wp = zeros(N, 1);

            case 'legendre'

                if debug; disp(' - Grid type: Legendre'); end

                J = diag(0.5 ./ sqrt(1 - (2 * (1:N-1)).^-2), 1) + diag(0.5 ./ sqrt(1 - (2 * (1:N-1)).^-2), -1);
                x = sort(eig(J));   % get Legendre points as eigenvalues of a Jacobi matrix

                x = (L(1) + L(2)) / 2 + (L(2) - L(1)) / 2 * x;
                w = ones(N, 1);
                wp = zeros(N, 1);

            case 'laguerre'

                if debug; disp(' - Grid type: Laguerre'); end

                J = diag(1:2:2*N-3) - diag(1:N-2, 1) - diag(1:N-2, -1);
                x = sort(eig(J));   % get Laguerre points as eigenvalues of a Jacobi matrix

                if isnan(decay)
                    decay = 2 * (L(2) - L(1)) / x(end);
                else
                    if debug; disp([' - Using decay scale D = ' num2str(decay)]); end
                end

                x = L(1) + [0; x * decay / 2];
                w = exp(- (x - L(1)) / decay);
                wp = -1 / decay * exp(-(x - L(1)) / decay);

            case 'hermite'

                if debug; disp(' - Grid type: Hermite'); end

                J = diag(sqrt((1:N-1) / 2), 1) + diag(sqrt((1:N-1) / 2), -1);
                x = sort(eig(J));   % get Hermite points as eigenvalues of a Jacobi matrix

                if isnan(decay)
                    decay = sqrt(2) * (L(2) - L(1)) / (2 * x(end));
                else
                    if debug; disp([' - Using decay scale D = ' num2str(decay)]); end
                end

                x = (L(1) + L(2)) / 2 + x * decay / sqrt(2);
                w = exp(- (x - (L(1) + L(2)) / 2).^2 / decay^2);
                wp = -2 / decay^2 * (x - (L(1) + L(2)) / 2) .* exp(- (x - (L(1) + L(2)) / 2).^2 / decay^2);

            case 'custom'

                if debug; disp(' - Grid type: Custom'); end

        end

        % Create  collocation matrix, M:
        D = x' - x;                                     % matrix of x_n - x_m
        D1 = D; D1(~D) = 1;                             % version of D with 0 -> 1 for performing products
        a_n = ones(N, 1 ) * (w' .* prod(D1));

        M = a_n' ./ (a_n .* D');
        M(1:(N+1):N^2) = sum(1 ./ D1) - 1 + (wp ./ w)'; % extra -1 removes the additional 1 gained from using D1

    otherwise

        error('Type:UnrecognisedType', ['The type ''' type ''' is not recognised.'])

end

% Flip x and M if flip = true:
if flip
    if debug; disp(' - Flipping segment'); end
    x = x(end:-1:1);
    M = M(end:-1:1,end:-1:1);
end

end