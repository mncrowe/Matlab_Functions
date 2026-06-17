function [s, x] = Gen_EVP(A, B, options)
% Solves the generalised EVP (s*A + B)*x = 0 for eigenvalue s and eigenvector x
%
% Inputs:
% - A,B: square matrices
% - options:
%   - n: number of eignevalues to find, default: 10
%   - s0: initial guess for eigenvalue, default: 0
%   - method: search method for eigs, default: 'lm'
%
% Note: this script finds the n eigenvalues closest to s0 by writing
% [(s0+s')*A + B]*x = 0  =>  [(1/s')*(s0*A+B) + A]*x = 0,
% and solving for eigenvalue 1/s'.

arguments
    A               (:,:)
    B               (:,:)
    options.n       (1,1) double {isinteger} = 10
    options.s0      (1,1) double             = 0
    options.method  (1,:) char               = 'lm'
end

% Determine eigenvalues and eigenvectors of modified system:

[V, S] = eigs(sparse(A), sparse(- B - options.s0 * A), options.n, options.method);

% Convert to eigenvalues of original system and rescale eigenvectors:

s = options.s0 + 1 ./ sum(S).';
x = reshape(V,[length(A) options.n]) ./ max(reshape(V, [length(A) options.n]));

% Remove NaN and inf and sort in descending order by real part:

x = x(:, and(~isnan(s), ~isinf(s)));
s = s(and(~isnan(s), ~isinf(s)));

[~, I] = sort(real(s), 'descend');
s = s(I);
x = x(:, I);

end

