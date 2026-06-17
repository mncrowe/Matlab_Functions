function n = sig_fig(n, f)
% Rounds the number n to f significant figures
%
% - n: number
% - f: number of figures

arguments
    n double
    f double {isinteger}
end

if n ~= 0
    s=n / abs(n);
    n=abs(n);
    n_exp = floor(log(n) / log(10));
    n = s * round(n * 10^(f-1) / 10^n_exp) * 10^n_exp / 10^(f-1);
end

end

