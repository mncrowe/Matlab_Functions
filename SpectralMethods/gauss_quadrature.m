function [I, x, w] = gauss_quadrature(f, options)
% Evaluates the integral of f using quadrature; int_D f(x) dx = sum w_i f(x_i)
%
% Inputs:
% - f: function or matrix. If a matrix, integration is performed through
%      the first dimension only and points should correspond to the values 
%      of the function at the underlying collication points
% - type: chacter array, the method used; options are 'chebyshev1', 
%         'chebyshev2' (default), 'legendre', 'laguerre', 'hermite', 
%         and 'custom'. Using 'custom' requires x and w to be specified
% - N: number of gridpoints used (default: 16)
% - L: integration domain (default: [-1 1] for hermite, [0 1] otherwise)
% - decay: decay rate parameter for 'laguerre' and 'hermite' methods
% - x, w: gridpoints and wieghts for custom quadrature scheme
%
% -------------------------------------------------------------------------
% Quadrature Types:
%
% Here p(x) is assumed to be a polynomial. All quadrature methods are exact
% when p(x) a polynomial of degree 2N-1. Let the integration domain be L =
% [a, b] and let c be the decay rate for Legendre and Hermite methods.
%
% Chebyshev-Gauss of the first kind ('chebyshev1'):
%
% This method used Chebyshev points of the first kind. It is appropriate
% for functions of the form f(x) = p(x) / sqrt((x-a)(b-x)).
%
% Chebyshev-Gauss of the second kind ('chebyshev2'):
%
% This method used Chebyshev points of the second kind. It is appropriate
% for functions of the form f(x) = p(x) * sqrt((x-a)(b-x)).
%
% Gauss-Legendre ('legendre'):
%
% This method uses Legendre points. It is appropriate for polynomials p(x).
%
% Gauss-Laguerre ('laguerre'):
%
% This method uses Laguerre points. It is appropriate for functions of the
% form f(x) = p(x) * exp(-cx). The decay rate, c, can be specified via
% 'decay' as decay = 1/c. If 'decay' is not specified, the decay rate will
% be set using the domain length, L, as c = x_N / [L(2)-L(1)].
% 
% Gauss-Hermite ('hermite'):
%
% This method uses Hermite points. It is appropriate for functions of the
% form f(x) = p(x) * exp(-(cx)^2). The decay rate, c, can be specified via
% 'decay' as decay = 1/c. If decay is not specified, the decay rate will be
% set using the domain length, L, as c = 2 * x_N / [L(2)-L(1)].
%
% Note: the Laguerre method implicitly assumes integration over [L(1) Inf]
% and the Hermite method implicitly assumes integration over [-Inf Inf]. L
% may be set to change the center of the domain or the decay scale via c.
% -------------------------------------------------------------------------

arguments
    f                   {functionormat}
    options.type  (1,:) char               = 'chebyshev2'
    options.N     (1,1) double {isinteger} = 16
    options.L     (1,:) double             = NaN
    options.decay (1,1) double             = NaN
    options.x     (:,1) double             = linspace(0, 1, 16)'
    options.w     (:,1) double             = ones(16, 1)
    options.flip  (1,1) logical            = false
end

% Get N as the length of f if f is a matrix:

if isa(f, 'double')
    N = size(f, 1);
else
    N = options.N;
end

% Set L based on size of input L value and method:

if numel(options.L) == 1

    if isnan(options.L)
        if strcmp(options.type, 'hermite')
            L = [-1 1];
        else
            L = [0 1];
        end
    else
        L = [-0.5 0.5] * options.L;
    end

else

    L = options.L;

end

% Define grid, x, and weights, w:

switch options.type

    case 'chebyshev1'

        x = (L(1) + L(2)) / 2 - (L(2) - L(1)) / 2 * cos((1/2:(N-1/2))*pi/N)';
        w = (L(2) - L(1)) / 2 * pi / N * sin((1/2:(N-1/2))*pi/N)';

    case 'chebyshev2'

        x = (L(1) + L(2)) / 2 - (L(2) - L(1)) / 2 * cos((0:(N-1))*pi/(N-1))';
        w = (L(2) - L(1)) / 2 * pi / (N - 1) * sin((0:(N-1))*pi/(N-1))';

    case 'legendre'

        J = diag(0.5 ./ sqrt(1 - (2 * (1:N-1)).^-2), 1) + diag(0.5 ./ sqrt(1 - (2 * (1:N-1)).^-2), -1);
        [V, D] = eig(J); % get Legendre points as eigenvalues of a Jacobi matrix
        [x, i] = sort(diag(D));
        x = (L(1) + L(2)) / 2 + (L(2) - L(1)) / 2 * x;
        w = (L(2) - L(1)) * (V(1, i).^2)';

    case 'laguerre'

        J = diag(1:2:2*N-3) - diag(1:N-2, 1) - diag(1:N-2, -1);
        [V, D] = eig(J); % get Laguerre points as eigenvalues of a Jacobi matrix
        [x, i] = sort(diag(D));

        if isnan(options.decay)
            decay = (L(2) - L(1)) / x(end);
        else
            decay = options.decay;
        end

        x = L(1) + [0; x * decay];
        w = 1 / decay * [0; (V(1, i).^2)'] .* exp((x - L(1)) / decay);

    case 'hermite'

        J = diag(sqrt((1:N-1) / 2), 1) + diag(sqrt((1:N-1) / 2), -1);
        [V, D] = eig(J); % get Laguerre points as eigenvalues of a Jacobi matrix
        [x, i] = sort(diag(D));

        if isnan(options.decay)
            decay = (L(2) - L(1)) / (2 * x(end));
        else
            decay = options.decay;
        end

        x = (L(1) + L(2)) / 2 + x * decay;
        w = sqrt(pi) * decay * ((V(1,i).^2)') .* exp((x - (L(1) + L(2)) / 2).^2 / decay^2);

    case 'custom'

        x = options.x;
        w = options.w;

    otherwise

        error('Type:UnrecognisedType', ['The type ''' type ''' is not recognised.'])

end

% Create vector f for function handle case:

if isa(f, 'function_handle')
    f = f(x);
end

% calculate integral

I = w' * f;

end

% Validation function for input f:

function functionormat(f)
    if ~isa(f, 'function_handle') && ~isa(f, 'double')
        error('Input:UnrecognisedType', 'The input f should be a function handle or matrix.')
    end
end

