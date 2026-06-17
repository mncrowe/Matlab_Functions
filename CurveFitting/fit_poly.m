function [coeff, func] = fit_poly(x, y, n)
% Fits an n^th degree polynomial to the given data using least squares
%
% - x: grid
% - y: data on grid
% - n: polynomial degree
%
% -------------------------------------------------------------------------
% Note: Polynomial interpolation is generally poor for high n though some
% choices of interpolation points, x, can work well, e.g. Chebyshev points.
% -------------------------------------------------------------------------

arguments
    x (:,1) double
    y (:,1) double
    n (1,1) double {mustBeInteger}
end

M = zeros(length(x), n + 1);

for in = 1:n+1
    M(:, in) = x .^ (in - 1);
end

coeff = M \ y;
func = @(z) (z .^ (0:n)) * coeff;

end

