function str_n = sig_fig_str(n,f)
% Rounds the number n to f significant figures and outputs as string
%
% - n: number
% - f: number of figures

arguments
    n double
    f int64
end

n = sig_fig(n,f);
str_n = num2str(n);

if n < 10^(f-1)
    if n == floor(n)
       str_n = [str_n '.'];
       for in = 1:(f - length(str_n) + 1)
           str_n = [str_n '0'];           
       end
    else
        if length(str_n) < f + 1
            for in = 1:(f - length(str_n) + 1)
                str_n = [str_n '0'];           
            end
        end
    end
end
    
end

