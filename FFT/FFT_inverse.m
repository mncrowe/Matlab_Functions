function f = FFT_inverse(f, d, Norm)
% Performs an inverse Fourier transform using the 'ifft' and 'fftshift'
% functions in the specified dimension(s).
%
% - f: array to transform (array)
% - d: dimensions to transform in (vector)
% - Norm: true - normalise (default), false - unscaled MATLAB tranform
%
% -------------------------------------------------------------------------
% Note: if Norm = true, the result is normalised by N where N is the
% product of the number of points in each dimension being tranformed in. 

arguments
    f          double
    d    (1,:) double {mustBeInteger} = 0
    Norm (1,1) logical                = true
end

N = size(f);
M = sum(N > 1);

if isequal(d, 0); d = 1:M; end

C = 1;

for i = d

    f = ifft(ifftshift(f, i), [], i);
    C = C * N(i);

end

if Norm; f = C * f; end

end

