function [x, k] = FFT_grid(N, L)
% Creates a (x,k) pair describing a grid in physical and Fourier space
%
% - N: number of gridpoints, integer
% - L: domain size;
%       scalar - x in [0 L]
%       vector - x in [L(1) L(2)]

arguments
    N (1,1) double {mustBeInteger}
    L (1,:) double
end

if numel(L) > 1
    X1 = L(1);
    X2 = L(2);
else
    X1 = 0; X2 = L;
end

dx = (X2 - X1) / N;
x = X1:dx:(X2-dx);
k = pi / dx * (-1:2/N:(1-2/N));

end

