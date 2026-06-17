function N = PeriodicConv2(M, C)
% Calculate the periodic convolution of M with C
%
% - M: input matrix, 2D array
% - C: convolution filter, 2D array (must be of Odd size in each dimension)

arguments
    M (:,:) double
    C (:,:) double {OddSize}
end

[Nx, Ny] = size(M);

c = size(C);

Cx = (c(1)-1)/2;
Cy = (c(2)-1)/2;

N = zeros(Nx, Ny);

for i = 1:Nx

    for j = 1:Ny

        Fx = mod((i-Cx:i+Cx)-1, Nx) + 1;
        Fy = mod((j-Cy:j+Cy)-1, Ny) + 1;

        N(i, j) = sum(C .* M(Fx, Fy), "all");

    end

end

end

function OddSize(C)
    [Cx, Cy] = size(C);
    if ~mod(Cx, 2) && ~mod(Cy, 2)
        eid = 'Size:NotOdd';
        msg = 'C must have odd length in each dimension.';
        error(eid, msg)
    end
end