function modon = CreateModon(grid, params)
% Creates a modon using the LCD or LRD solution

arguments
    grid              struct
    params.U    (1,1) double              = 1
    params.a    (1,1) double              = 1
    params.beta (1,1) double              = 0
    params.R    (1,1) double              = inf
    params.mode (1,1) double {isinteger}  = 1
end

% Get polar coordinates:

r = sqrt(grid.x.^2 + grid.y.^2);
theta = atan2(grid.y, grid.x);

r(r == 0) = eps;    % avoid inf values when dividing by r

% Define intermediate parameters and functions:

U = params.U;
a = params.a;
beta = params.beta;
R = params.R;

p = sqrt(beta/U + 1/R^2);

J1 = @(x) besselj(1, x);
J1p = @(x) (besselj(0, x) - besselj(2, x)) / 2;
K1 = @(x) besselk(1, x);
K1p = @(x) (-besselk(0, x) - besselk(2, x)) / 2;

% Find the required root of the Bessel function J1(x):

z = (params.mode + 1/4) * pi;
K0 = z - 3 / (8 * z) + 12 / (8 * z)^3;  % use McMahon's expansion to estimate root
K0 = fzero(J1, K0);                     % refine estimate using root finding

% Define radial structure functions:

if p == 0

    % Create LCD:

    type = 'LCD';

    % Set K using first root of Bessel function J_1:

    K = K0 / a;

    % Define coefficients for psi and q:

    A = -U * a^2;
    B = 2 * U / (K * J1p(K * a));

    % Define radial structure functions for psi and q:

    Psi   = @(r) A ./ r .* (r >= a) + (B * J1(K * r) - U * r) .* (r < a);
    Psi_r = @(r) -A ./ r.^2 .* (r >= a) + (B * K * J1p(K * r) - U) .* (r < a);
    Q     = @(r) -K^2 * B * J1(K * r) .* (r < a);
    Q_r   = @(r) -K^3 * B * J1p(K * r) .* (r < a);

else

    % Create LRD:

    type = 'LRD';

    % Find K using root finding:

    f = @(x) x .* J1p(x) - (1 + x.^2 / (p^2 * a^2)) .* J1(x) + x.^2 .* J1(x) * K1p(p * a) / (p * a * K1(p * a));
    Kp = fzero(f, K0);      % find Kp by looking for roots near K0, mode may differ for large beta or small R
    K = a * sqrt(Kp^2 + 1/R^2);

    % Define coefficients for psi and q:

    A = -U * a / K1(p * a);
    B = p^2 * U * a / (Kp^2 * J1(Kp * a));

    % Define radial structure functions for psi and q:

    Psi   = @(r) A * K1(p * r) .* (r >= a) + (B * J1(Kp * r) - U * (Kp^2 + p^2) / Kp^2 * r) .* (r < a);
    Psi_r = @(r) A * p * K1p(p * r) .* (r >= a) + (B * Kp * J1p(Kp * r) - U * (Kp^2 + p^2) / Kp^2) .* (r < a);
    Q     = @(r) beta / U * Psi(r) .* (r >= a) - (K^2 / a^2 * Psi(r) + (U * K^2 / a^2 + beta) * r) .* (r < a);
    Q_r   = @(r) beta / U * Psi_r(r) .* (r >= a) - (K^2 / a^2 * Psi_r(r) + (U * K^2 / a^2 + beta)) .* (r < a);

end

% Define psi, q, psi_x, psi_y, q_x, q_y:

psi   = Psi(r) .* sin(theta);
psi_x = (Psi_r(r) - Psi(r) ./ r) .* sin(theta) .* cos(theta);
psi_y = Psi_r(r) .* sin(theta).^2 + Psi(r) ./ r .* cos(theta).^2;
q     = Q(r) .* sin(theta);
q_x   = (Q_r(r) - Q(r) ./ r) .* sin(theta) .* cos(theta);
q_y   = Q_r(r) .* sin(theta).^2 + Q(r) ./ r .* cos(theta).^2;

% Assign fields to output structure:

if isa(grid.K2, 'gpuArray')
    modon.psi = gpuArray(psi);
    modon.q = gpuArray(q);
    modon.psi_x = gpuArray(psi_x);
    modon.psi_y = gpuArray(psi_y);
    modon.q_x = gpuArray(q_x);
    modon.q_y = gpuArray(q_y);
else
    modon.psi = psi;
    modon.q = q;
    modon.psi_x = psi_x;
    modon.psi_y = psi_y;
    modon.q_x = q_x;
    modon.q_y = q_y;
end

modon.K = K;
modon.U = U;
modon.a = a;
modon.beta = beta;
modon.R = R;
modon.type = type;

end