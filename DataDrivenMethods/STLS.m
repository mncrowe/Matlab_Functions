function x = STLS(A, b, l, N)
% Implements Sequentially Thresholded Least Squares
% Finds sparse x which minimises Ax = b, approximating the loss function 
% L(x) = |b-Ax|_2^2 + l*|x|_0
%
% Inputs:
% - A: matrix, n x m
% - b: matrix, n x d
% - l: lambda (sparsification parameter), default: 0.1
%   - scalar; la - absolute threshold parameter
%   - vector; [la lr], lr - relative threshold parameter
% - N: maximum number of iterations, default: 10
%
% Outputs:
% - x: matrix, m x d
%
% -------------------------------------------------------------------------
% Notes: the code will terminate if no additional elements are set to zero
% by the thresholding step. Based on code by S. L. Brunton.
%
% This function is intended for use with SINDy. The SINDy coefficients are
% determined as (for example) c = STLS(Theta, dX_dt, 0.2) where Theta is
% the library matrix in which the library elements appear as columns and
% dX_dt is a matrix with the time-derivatives as columns.

arguments
    A (:,:) double
    b (:,:) double
    l (1,1) double                 = 0.1
    N (1,1) double {mustBeInteger} = 10
end

% Set relative threshold to zero if not given:
if numel(l) == 1; l = [l 0]; end

% Check dimensions are consistent and define n:
n = CheckDims(A, b);

% Initial guess using linear regression:
x = A \ b;
spel = false(size(x));

% Iterative regression, set parameters small than l to 0 in each loop:
for iN = 1:N
    if isequal(spel,(abs(x) < (l(1) + l(2) * max(abs(x), [], 'all'))))
        break;  % exit if number of sparse elements has not changed
    else
        spel = (abs(x) < (l(1) + l(2) * max(abs(x), [], 'all'))); % elements to zero
    end
    x(spel) = 0;
    for in = 1:n
        x(~spel(:, in), in) = A(:, ~spel(:, in)) \ b(:, in);
    end
end

end

function n = CheckDims(A, b)
    [A1, ~] = size(A);
    [B1, B2] = size(b);
    if A1 == B1
        n = B2;
    else
        eid = 'Size:DimensionMismatch';
        msg = 'A and b must have the same length in the first dimension.';
        error(eid, msg)
    end
end