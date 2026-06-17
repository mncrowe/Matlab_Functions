function F_avg = window_average(F, t, dim, options)
% Time averages a given input over a given window
%
% - F: field to average
% - t: time vector (default: 1:N)
% - dim: dimension to average over (default: length(size(F)), final dimension)
% - options:
%   - range: vector [i1 i2], indices of time points where average is calculated (default: [1 Nt])
%   - window: vector [w1 w2], window to average over, relative to each index (default: [-1 1])
%   - method: averaging method;
%       - 0: central Riemann sum/trapezoidal rule (default)
%       - 1: left Riemann sum
%       - 2: right Riemann sum

arguments
    F                    double
    t              (1,:) double                     = 1:size(F, length(size(F)))
    dim            (1,1) double {CheckDim(dim,F,t)} = length(size(F))
    options.range  (1,:) double                     = [1 length(t)]
    options.window (1,:) double {mustBeInteger}     = [-1 1]
    options.method (1,1) double {mustBeInteger}     = 0
end

% Get size of F and permute so averaged dimension is last:

N = length(size(F));        % number of dimensions

I = 1:N; I(dim) = [];       % index of non-differentiated dimensions
F = permute(F, [I dim]);

S = size(F);
Nt = S(end);

F = reshape(F, [], Nt);

% Get output size from averaging range and define output: 

Nr = options.range(2) - options.range(1) + 1;

F_avg = zeros(prod(S(1:end-1)), Nr);

% Loop through averaging dimension performing average:

for i = options.range(1):options.range(2)

    % Define indices corresponding to edges of averaging window:

    j1 = max(1, i + options.window(1));
    j2 = min(Nt, i + options.window(2));

    % Define dt over averaging window:

    dt = reshape(t(j1+1:j2) - t(j1:j2-1), [1 j2-j1]);
    T = t(j2) - t(j1);

    % Define and apply averaging filter:

    switch options.method

        case 0 % trapezoidal rule

            M = ([dt 0] + [0 dt]) / 2;

        case 1 % left Riemann sum

            M = [dt 0];

        case 2 % right Riemann sum

            M = [0 dt];

        otherwise % not implemented error

        eid = 'Method:NotImplemented';
        msg = ['Using method = ' num2str(options.method) ' has not been implemented.'];
        error(eid, msg)

    end

    F_avg(:, i - options.range(1) + 1) = sum(F(:, j1:j2) .* M, 2) / T;

end

% Reshape and permute output to match required size and dimension order:

F_avg = reshape(F_avg, [S(1:end-1) Nr]);
F_avg = permute(F_avg, [I(I<dim) N I(I>dim)]);

end

function CheckDim(dim, F, t)
    mustBeInteger(dim)
    if length(t) ~= size(F, dim)
        eid = 'TimeVector:InconsistentLength';
        msg = ['The length of t does not match the length of dimension ' num2str(dim) ' of F.'];
        error(eid, msg)
    end

end