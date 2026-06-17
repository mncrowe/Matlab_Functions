function sol = timestep(options)
% Solve the system: D*d(phi)/dt = L*phi + N(phi,t) using implicit-explicit (IMEX) Runge-Kutta (RK) schemes
%
% - D: LHS matrix, n x n
% - L: linear RHS matrix, n x n
% - N: nonlinear RHS, anonymous function, N = @(phi,t) ..., enter constant forcing as N = @(phi,t) [v1; ...; vn]
% - dt: timestep, assumed to be constant, (set to min(dt, diff(t)))
% - Nt: number of timesteps, if t is specified, this will be overwritten using Nt = ceil(t(end) / dt)
% - Ns: number of evenly spaced time saves, if t is specified, this will be overwritten
% - t: vector of time saves, if not specified, t = linspace(0, Nt*dt, Ns+1)
% - phi0: initial condition for phi
% - method: timestepping method used;
%       'TRAP'   - 1 stage, IMEX trapezoidal method, 1st order DIRK+ERK scheme (default)
%       'RK111'  - 1 stage, RK111, 1st order DIRK+ERK scheme (implicit backward/explicit forward Euler)
%       'RK222'  - 2 stage, RK222, 2nd order DIRK+ERK scheme (Ascher et. al. 97)
%       'RK332'  - 3 stage, RK332, 2nd order DIRK+ERK scheme (Koto 06)
%       'RKSMR'  - 3 stage, RKSMR, 3-eps order DIRK+ERK scheme (Spalart et. al. 91)
%       'RK443'  - 4 stage, RK443, 3rd order DIRK+ERK scheme (Ascher et. al. 97)
%       'EULER'  - 1 stage, RK1, 1st order ERK (forward Euler) (not recommended for singular D)
%       'RK2'    - 2 stage, RK2, 2nd order ERK scheme (not reccommended for singular D), method depends on alpha
%       'RK3'    - 3 stage, RK3, 3rd order ERK scheme (not reccommended for singular D), method depends on alpha and beta
%       'SSPRK3' - 3 stage, Strong Stability Preserving RK3, 3rd order ERK scheme
%       'RK4'    - 4 stage, RK4, 4th order ERK scheme (not reccommended for singular D)
%       'CUSTOM' - custom, enter cell array, C, with C{1} = s, C{2} = H, C{3} = A, C{4} = c
% - alpha: alpha value for arbitrary RK2 and RK3 methods (default: 1/2)
% - beta: beta value for arbitrary RK3 methods (default: 1)
% - C: custom method; enter cell array of C{1:4} = {s, H, A, c} for custom RK scheme, required for method = 'RK2'
% - debug: if true, progress information is printed to screen (default: true)
% - savetol: tolerance used to determine if current time is close enough to save time when saving sol
% - info: if true, additional solver information is included in sol
% - exact_t: if true, sol.t is overwritten with the value of t corresponding to sol.phi, may lead to unevenly spaced sol.t
%
% ----------------------------------------------------------------------------
% Note: Terminology; DIRK (diagonally implicit Runge-Kutta), ERK (explicit
%       Runge-Kutta), IMEX (implicit-explicit), SMR (Spalart, Moser & Rogers)
% 
% Note: alpha and beta can be specified to choose different RK2 and RK3
%       methods. For RK2: alpha = 1 (Heun's method), 1/2 (midpoint method, 
%       default), 2/3 (Ralston method). For RK3: (alpha, beta) = (1/2, 1)
%       (Kutta's method, default), (1/3, 2/3) (Heun's method), (1/2, 3/4)
%       (Ralston's method), (8/15, 2/3) (Wray's method). Note SSPRK3 has
%       (alpha, beta) = (1, 1/2) but is included seperately for ease of
%       use. The values of alpha and beta are saved to sol if info = true.
%
% Note: H, A, c do not use the full Butcher tableau as sometimes presented.
%       A is the s x s bottom left square block, H is the 2nd to final rows and
%       c is the first s entries only (0 to s-1), so c_0 = 0. This is as any
%       remaining values are zeros or ones and are known. See examples below.
%       Therefore: sum_j A_ij = sum_j H_ij = c_i, for i=1..s.
%       Here c_s = 1 is not included in c. A and H must be lower triangular.
%
% Note: if the required save points, t, do not line up with the calculation
%       time (multiples of dt) then the solution at time t will be saved to
%       sol.phi(:, i) when t+savetol*dt > t.sol(i). Setting exact_t = true
%       will update sol.t(i) with the actual t value.
% ----------------------------------------------------------------------------

arguments
    options.D       (:,:) double             = 1
    options.L       (:,:) double             = 1
    options.N             function_handle    = @(psi, t) 0;
    options.dt      (1,1) double             = 0.1
    options.Nt      (1,1) double {isinteger} = 100
    options.Ns      (1,1) double {isinteger} = 100
    options.t       (1,:) double             = 0
    options.phi0    (:,1) double             = 1
    options.method  (1,:) char               = 'TRAP'
    options.alpha   (1,1) double             = 0.5
    options.beta    (1,1) double             = 1
    options.C             struct
    options.debug   (1,1) logical            = true
    options.savetol (1,1) double             = 1e-2
    options.info    (1,1) logical            = true
    options.exact_t (1,1) logical            = false
end

% Define timestepping parameters:

if numel(options.t) > 1
    sol.t = options.t;
    dt = min(options.dt, min(diff(options.t)));
    Nt = ceil((options.t(end) - options.t(1)) / dt);
else
    sol.t = linspace(0, options.Nt * options.dt, options.Ns + 1);
    dt = options.dt;
    Nt = options.Nt;
end

% Get problem size:

NL = max([length(options.L), length(options.D), length(options.phi0)]);

% Define IMEX Runge-Kutta tableaus, H is s x (s+1), A is s x s, c is 1 x s:

a = options.alpha;
b = options.beta;

switch options.method

    case 'TRAP'

        s = 1;
        H = [0.5 0.5];
        A = 1;
        c = 0;

    case 'RK111'

        s = 1;
        H = [0 1];
        A = 1;
        c = 0;

    case 'RK222'

        s = 2;
        H = [0 1-1/sqrt(2) 0; 0 1/sqrt(2) 1-1/sqrt(2)];
        A = [1-1/sqrt(2) 0; -1/sqrt(2) 1+1/sqrt(2)];
        c = [0 1-1/sqrt(2)];

    case 'RK332'

        s = 3;
        H = [0 1 0 0; 0 -1/2 1 0; 0 -1 1 1];
        A = [1 0 0; 1/2 0 0; 0 0 1];
        c = [0 1 1/2];

    case 'RKSMR'

        s = 3;
        H = [29/96 37/160 0 0; 29/96 5/32 5/24 0; 29/96 5/32 3/8 1/6];
        A = [8/15 0 0; 1/4 5/12 0; 1/4 0 3/4];
        c = [0 8/15 2/3];

    case 'RK443'

        s = 4;
        H = [0 1/2 0 0 0; 0 1/6 1/2 0 0; 0 -1/2 1/2 1/2 0; 0 3/2 -3/2 1/2 1/2];
        A = [1/2 0 0 0; 11/18 1/18 0 0; 5/6 -5/6 1/2 0; 1/4 7/4 3/4 -7/4];
        c = [0 1/2 2/3 1/2];

    case 'EULER'

        s = 1;
        H = [1 0];
        A = 1;
        c = 0;

    case 'RK2'

        s = 2;
        H = [a 0 0; 1-1/(2*a) 1/(2*a) 0];
        A = [a 0; 1-1/(2*a) 1/(2*a)];
        c = [0 a];

    case 'RK3'

        s = 3;
        H = [a 0 0 0; b/a*(b-3*a*(1-a))/(3*a-2) -b/a*(b-a)/(3*a-2) 0 0; 1-(3*a+3*b-2)/(6*a*b) (3*b-2)/(6*a*(b-a)) (2-3*a)/(6*b*(b-a)) 0];
        A = [a 0 0; b/a*(b-3*a*(1-a))/(3*a-2) -b/a*(b-a)/(3*a-2) 0; 1-(3*a+3*b-2)/(6*a*b) (3*b-2)/(6*a*(b-a)) (2-3*a)/(6*b*(b-a))];
        c = [0 a b];

    case 'SSPRK3'

        s = 3;
        H = [1 0 0 0; 1/4 1/4 0 0; 1/6 1/6 2/3 0];
        A = [1 0 0; 1/4 1/4 0; 1/6 1/6 2/3];
        c = [0 1 1/2];

    case 'RK4'

        s = 4;
        H = [1/2 0 0 0 0; 0 1/2 0 0 0; 0 0 1 0 0; 1/6 1/3 1/3 1/6 0];
        A = [1/2 0 0 0; 0 1/2 0 0; 0 0 1 0; 1/6 1/3 1/3 1/6];
        c = [0 1/2 1/2 1];

    case 'CUSTOM'

        s = options.C{1};
        H = options.C{2};
        A = options.C{3};
        c = options.C{4};

    otherwise
        error('Method:UnrecognisedValue', 'The method entered is not recognised.')
end

% Decompose LHS implicit terms for fast inversion, convert D and L to sparse:

warning('off','MATLAB:nearlySingularMatrix')    % large matrices usually are ill-conditioned by MATLAB's definition, instead we check for NaNs

if options.debug; disp('Decomposing matrices...'); end

for i = 1:s
    if options.debug; disp([num2str(i) ' of ' num2str(s)]); end
    M{i} = decomposition(options.D - dt * H(i,i+1) * options.L);
end

D = sparse(options.D);
L = sparse(options.L);

% Define solution, define intermediate stages, K, and set initial phi and time:

sol.phi = zeros(NL, length(sol.t));
sol.phi(:, 1) = options.phi0;

K = zeros(NL, s+1);

phi = options.phi0;  % set current solution value
t = sol.t(1);        % set current time

i_save = 2;          % position in dim 2 of sol.phi where next saved field save should go

% Timestep system using RK method:

if options.debug; disp('Starting timestepping...'); end

runtime = tic;   % start timer

for it = 1:Nt

    % IMEX Runge-Kutta loop:
    
    K(:, 1) = phi;

    for i = 1:s                                         % loop through RK stages
        RHS = zeros(NL, 1);                             % (re)set RHS to zero
        for j = 1:i                                     % build RHS terms
            RHS = RHS + H(i, j) * L * K(:, j) + A(i, j) * options.N(K(:, j), t + c(j) * dt);
        end
        K(:, i+1) = M{i} \ (D * K(:, 1) + dt * RHS);    % solve implicit system for next K 
    end

    % Update solution values:

    phi = K(:, s+1);                                    % update phi_n -> phi_{n+1}
    t = t + dt;                                         % update t -> t + dt

    % Save solution if current time is close to a time-point in the output t:

    if i_save <= length(sol.t) && (t + dt * options.savetol) >= sol.t(i_save)

        sol.phi(:, i_save) = phi;
        if options.exact_t; sol.t(i_save) = t; end

        i_save = i_save + 1;

        % Print progress:

        if options.debug
            counter(toc(runtime), it / Nt, 0, ['t = ' num2str(t)])
        end

        % Display error if NaN or Inf in solution:

        if max(isnan(phi), isinf(phi))
            error('Value:NaN', ['NaN error, t = ' num2str(t)]);
        end
            
    end

end

% Print total runtime:

if options.debug
    disp(['Timestepping complete, elapsed time = ' num2str(toc(runtime))])
end

% Save additional solver information to sol:

if options.info
    
    sol.runtime = toc(runtime);
    sol.method = options.method;
    sol.dt = dt;
    
    if strcmp(options.method, 'RK2')
        sol.alpha = a;
    end

    if strcmp(options.method, 'RK3')
        sol.alpha = a;
        sol.beta = b;
    end
    
    if strcmp(options.method, 'CUSTOM')
        sol.C = options.C;
    end

end

end