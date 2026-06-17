function [Fmn, a, b, expQ, expR] = Pade_2D(F, x, y, m, n, round)
% Calculates the approximant P = Q(x,y)/R(x,y) where Q is order (m1, m2) in
% (x, y) and R is order (n1, n2) in (x, y)
%
% Inputs:
% - f: function to approximate, vector
% - x: position coordinate in dimension 1, vector or matrix (default: 1:Nx)
% - y: position coordinate in dimension 2, vector or matrix (default: 1:Ny)
% - m: vector, polynomial orders of numerator (default: [2, 2])
% - n: vector, polynomial orders of denominator (default: [2, 2])
% - round: round coefficients to given tolerance (default: 0)
%
% Outputs:
% - Fmn: anonymous function representing Pade approximant
% - a: numerator coefficients
% - b: denominator coefficients
% - expQ, expR: exponents of each term in Q and R, column are [x y] indices
%
% -------------------------------------------------------------------------
% Notes:
% Setting round = 0 performs no rounding of the coefficients, a value of
% round > 0, rounds the coefficients to that number of significant figures.
%
% x and y may be entered as 2D matrices with varying x corresponding to 
% dimension 1 of each and varying y to dimension 2. The values of F should 
% correspond to the values at these points hence x, y, F must all be the 
% same size if matrix inputs are used.
%
% Enter m and/or n as vectors of length 2 or 3. If entered as length 3, the
% final element gives the maximum order of each term. For example m = [2,2]
% gives Q = a + bx + cx^2 + dy + exy + fx^2y + gy^2 + hxy^2 + gx^2y^2 while
% m = [2,2,2] gives Q = a + bx + cx^2 + dy + exy + fy^2.
%
% This function requires 'sig_fig.m'.
% -------------------------------------------------------------------------

arguments
    F     (:,:) double
    x     (:,:) double {CheckSize(1,x,F)} = 1:size(F,1);
    y     (:,:) double {CheckSize(2,y,F)} = 1:size(F,2);
    m     (1,:) double {mustBeInteger}    = [2 2]
    n     (1,:) double {mustBeInteger}    = [2 2]
    round (1,1) double {mustBeInteger}    = 0
end

% Get (Nx, Ny) from size of F:

Nx = size(F, 1);
Ny = size(F, 2);

% Reshape x and y to 2D arrays if required:

if numel(x) == Nx
    x = reshape(x, [Nx 1]) * ones(1, Ny);
end

if numel(y) == Ny
    y = ones(Nx, 1) * reshape(y, [1, Ny]);
end

% Define arrays of exponents in x and y for polynomials Q and R:

expQ = zeros(m(1) + 1, m(2) + 1, 2);
expQ(:, :, 1) = (0:m(1))' * ones(1, m(2) + 1);
expQ(:, :, 2) = ones(m(1) + 1, 1) * (0:m(2));

expR = zeros(n(1) + 1, n(2) + 1, 2);
expR(:, :, 1) = (0:n(1))' * ones(1, n(2) + 1);
expR(:, :, 2) = ones(n(1) + 1, 1) * (0:n(2));

% Reshape exponent arrays to two vectors representing x and y powers:

expQ = reshape(expQ, [], 2);
expR = reshape(expR, [], 2);

% Remove terms that exceed maximum total order:

if length(m) == 2
    m(3) = m(1) + m(2);
end

if length(n) == 2
    n(3) = n(1) + n(2);
end

expQ(expQ(:, 1) + expQ(:, 2) > m(3), :) = [];
expR(expR(:, 1) + expR(:, 2) > n(3), :) = [];
expR(1, :) = []; % R has 1 as it's constant term by assumption

% Reshape x, y and F arrays to vectors:

xp = reshape(x, [Nx*Ny 1]);
yp = reshape(y, [Nx*Ny 1]);
Fp = reshape(F, [Nx*Ny 1]);

% Calculate matrix of x^n * y^m terms for linear regression step:

M = [xp.^(expQ(:, 1)').*yp.^(expQ(:, 2)') -Fp.*xp.^(expR(:, 1)').*yp.^(expR(:, 2)')];

% Perform linear regression step:

C = M \ Fp;

% Round coefficients if required:

if round > 0
    for i = 1:length(C)
        C(i) = sig_fig(C(i), round);
    end
end

% Define output coefficient arrays:

a = C(1:length(expQ));
b = [1; C((length(expQ) + 1):(length(expQ) + length(expR)))];
expR = [0 0; expR];

% Define output Pade approximant function:

Fmn = @(x,y) reshape((reshape(x, [Nx*Ny 1]) .^ (expQ(:, 1)') .* reshape(y, [Nx*Ny 1]) .^ (expQ(:, 2)')) * a ./ ...
    ((reshape(x, [Nx*Ny 1]) .^ (expR(:, 1)') .* reshape(y, [Nx*Ny 1]) .^ (expR(:, 2)')) * b), [Nx Ny]);

end

function CheckSize(dim, x, F)
    if ~isequal(size(x), size(F)) && ~(isequal(length(x), size(F, dim)) & numel(x) == length(x))
        eid = 'Length:NotEqual';
        msg = 'Vectors x and f must have the same length.';
        error(eid, msg)
    end

end