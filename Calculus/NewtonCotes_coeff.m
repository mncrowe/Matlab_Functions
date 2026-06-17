function H = NewtonCotes_coeff(n)
% Calculates the coefficients for an n-point Newton-Cotes rule
%
% - n: number of points

arguments
    n (1,1) double {mustBeInteger} = 2
end

H = zeros(1, n);

syms t r

P = symprod(t - r, r, 0, n - 1);

for i = 0:n-1
    I = int(P / (t-i), t, 0, n - 1);
    H(i + 1) = (-1)^(n - i - 1) / (factorial(i) * factorial(n-i-1)) * I;
end

end