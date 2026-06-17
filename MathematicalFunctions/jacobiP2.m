function P = jacobiP2(n, a, b, x)
% A better Jacobi polynomial J_n^{a,b}(x) that's less slow than the Matlab
% one though only works for integer a and b
%
% - n, a, b: degree and parameters for polynomial, integers
% - x: grid points where polynomial is evaluated, array

arguments
    n (1,1) {mustBeInteger}
    a (1,1) double
    b (1,1) double
    x double
end

P = 0;

for  s = 0:n
    P = P + 1 / (GammaFact(s) * GammaFact(n+a-s) * GammaFact(b+s) * GammaFact(n-s)) ...
        * ((x-1) / 2).^(n-s) .* ((x+1) / 2).^s;
end

P = P * GammaFact(n+a) * GammaFact(n+b);

end

function a = GammaFact(a)
% Calculates a! using factorial if a is an integer and the Gamma function
% if a non-integer
%
% - a: input number

if mod(a, 1) == 0
    a = factorial(a);
else
    a = gamma(a + 1);
end

end
