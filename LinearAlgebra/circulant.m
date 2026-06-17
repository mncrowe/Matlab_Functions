function M = circulant(N, i)
% Creates a circulant matrix consisting of a single (looped) diagonal with entry 1
%
% - N: Matrix size (square)
% - i: Index of diagonal, e.g. 0; leading diagonal, 1; diagonal above 0, -1; diagonal below 0, etc

arguments
    N (1,1) int64
    i (1,1) int64
end

i = rem(i, N);
l = N - abs(i);

if i == 0
    M = eye(N);
else
    M = diag(ones(1, l), i) + diag(ones(1, abs(i)), -sign(i) * l);
end

end

