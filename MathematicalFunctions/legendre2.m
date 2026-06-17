function out = legendre2(m, n, x)
% Outputs the values of the associated legendre polynomial, P^m_n, at point(s) x
%
% - m: order of polynomial, integer, may be negative
% - n: degree of polynomial, integer, may be negative
% - x: vector of points

arguments
    m (1,1) {mustBeInteger}
    n (1,1) {mustBeInteger}
    x double
end

A = 1;

if n < 0; n= -n-1; end
if m < 0; m = -m; A = (-1)^m * factorial(n-m) / factorial(n+m); end

leg = legendre(n, x);
leg = reshape(leg, n+1, numel(leg) / (n+1));

out = A * reshape(squeeze(leg(m+1, :)), size(x));

end