function [x, n] = newton_raphson(f, x_0, options)
% solves f(x) = 0 with initial guess x_0 for f:R^n -> R^n using a
% generalised Newton-Raphson method
%
% Inputs:
% - f: function, e.g. f(x) = @(x) [x(1).^2+x(2).^2; x(1)]
% - x_0: initial guess of root, vector
% - eps_0: error tolerance (default: 1e-3)
% - N: max number of iterations (default: 100)
% - delta: grid spacing for numerical differentiation (default: 1e-2)
% - output: true - show f(x) at each step, false - don't show (default: true)
%
% Outputs:
% - x: root
% - n: number of iterations

arguments
    f              function_handle
    x_0            (:,1) double
    options.eps0   (1,1) double          = 1e-3
    options.N      (1,1) {mustBeInteger} = 100
    options.delta  (1,1) double          = 1e-2
    options.output (1,1) logical         = true
end

dim = length(x_0);
x = x_0;
n = 0;

while f(x)' * f(x) > options.eps0^2
    
    M = zeros(dim);
    for i = 1:dim
        delta_vec = options.delta * double(1:dim == i)';
        M(:, i) = (f(x+delta_vec) - f(x-delta_vec)) / (2 * options.delta);
    end
    x = x - M^(-1) * f(x);
    n = n+1;
    if n == options.N
        error(['No convergence in ' num2str(options.N) ' iterations'])
    end
    
    if options.output
        disp(['Iteration ' num2str(n) ':'])
        disp([' - |f(x)|^2 = ' num2str(f(x)' * f(x))])
    end

end

end