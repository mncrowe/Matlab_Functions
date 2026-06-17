function [p, c] = fit_polynomial(x, f, m)
% Fits the function f(x) to an order m polynomial, p(x)
%
% Inputs:
% - x: grid: vector of length n
% - f: function values, function handle or vector of length n
% - m: order of polynomial, (default: n-1)
%
% Outputs:
% - p: inline function, polynomial of degree m, g = sum_{i=0}^m [g_i x^i]
% - c: coefficients of p, [g_0, g_1, ..., g_{n-1}, g_n]

arguments
    x (:,1) double
    f       {doubleorfunction(x,f)}
    m (1,1) double {isinteger}      = length(x)-1
end

if isa(f, 'function_handle')
    f = f(x);
else
    f = reshape(f, [], 1);
end

% Solve for coefficients using least squares method:

M = (x * ones(1, m + 1)) .^ (0:m);
c = M \ f;
c(abs(c) < 1e2 * eps * max(abs(c))) = 0;    % set small entries to zero

% Create p(x) using symbolic maths:

syms z
p = sum(c' .* z .^(0:m));
p = matlabFunction(p);

end

% Validation for function or double input:

function doubleorfunction(f, x)
    if ~(isa(f, 'double') && length(f) == length(x)) && ~isa(f, 'function_handle')
        error('TypeError:UnexpectedType', 'The function f should be a vector of the same length as x or a function handle.')
    end
end

