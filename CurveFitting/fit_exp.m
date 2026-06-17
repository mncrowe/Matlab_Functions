function [a, b] = fit_exp(x, y)
% Fits the data to a curve y = a*exp(b*x) and returns [a, b]
%
% - x: independent variable
% - y: dependent variable

arguments
    x (1,:) double
    y (1,:) double
end

c = polyfit(x, log(y), 1);

a = exp(c(2));
b = c(1);

end

