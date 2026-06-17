function [P, C, I, R] = Orthog_Poly(N, a, b, w, err, norm)
% Finds the first N + 1 orthogonal polynomials on [a, b] w.r.t. the weight
% function w(x) where <P_n, P_m> = integral(P_n(x) P_m(x) w(x), a, b)
%
% Inputs:
% - N: order of highest polynomial to calculate (integer), default: 1
% - [a, b]: domain to integrate over, may be infinite, default: [0 1]
% - w: weight function w(x), function, default: @(x) 1
% - err: error in numerical integration, default: 1e-6
% - norm: normalisation method, "left": P_n(a) = 1, "right": P_n(b) = 1, "int" <P_n, P_n> = 1
%
% Outputs:
% - P: orthogonal polynomial as a function P(n, x), accepts 2D array inputs
% - C: matrix of coefficients, C(i, :) is the cofficients of P_i
% - I: matrix of integral values of x^n w(x)
% - R: matrix containing the roots of each P_n(x) as rows
%
% Note: roots are found by finding the eigenvalues of the companion matrix.
% This method is not optimal for finding roots of some sets of orthogonal 
% polynomials where some coefficients are particularly large or small, e.g.
% Laguerre polynomials. In such cases, other methods (such as finding the 
% eigenvalues of a corresponding Jacobi matrix) may be better.

arguments
    N       (1,1)   double {isinteger}  = 1
    a       (1,1)   double              = 0
    b       (1,1)   double              = 1
    w       (1,1)   function_handle     = @(x) 1
    err     (1,1)   double              = 1e-6
    norm    (1,1)   string              = "left"    % left, right, int
end

% build Hankel matrix I by integrating powers of x w.r.t. w(x):

I = zeros(N + 1);

for i = 1:2*N+1
    I = I + diag(ones(min(i, 2*N + 2 - i), 1), i - N - 1) * integral(@(x) x.^(i-1) .* w(x), a, b, AbsTol = err, RelTol = err);
end

I = I(end:-1:1, :);

% create coefficient matrix, C:

C = eye(N + 1);     % set P_0(x) = 1

for n = 1:N

    % build linear system for coefficients of P_n given P_i for all i < n:
    
    M = C(1:n, 1:n) * I(1:n, 1:n);
    c = -C(1:n, 1:n) * I(1:n, n+1);

    % solve for coefficients of next polynomial:

    C(n + 1, 1:n) = M \ c;

    switch norm
        case "left"
            C(n + 1, :) = C(n + 1, :) / (a.^(0:n) * C(n + 1, 1:n+1)');
        case "right"
            C(n + 1, :) = C(n + 1, :) / (b.^(0:n) * C(n + 1, 1:n+1)');
        case "int"
            C(n + 1, :) = C(n + 1, :) / sqrt(C(n + 1, 1:n+1) * I(1:n+1, 1:n+1) * C(n + 1, 1:n+1)');
        otherwise
            error("value of 'norm' is invalid.")
    end

end

% define orthogonal polynomials to accept 2D array inputs:

P = @(n, x) sum(reshape(C(n + 1, 1:n+1), 1, 1, []) .* x .^ reshape(0:n, 1, 1, []), 3);

% find the roots of the n^th polynomial:

if nargout > 3

    R = zeros(N);

    for i = 1:N

        A = diag(ones(i - 1, 1), -1);
        A(1,:) = -C(i + 1, i:-1:1) / C(i + 1, i + 1);
        R(i, 1:i) = sort(eig(A));

    end

end

end