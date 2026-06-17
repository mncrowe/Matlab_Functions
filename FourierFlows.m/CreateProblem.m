function prob = CreateProblem(options)

arguments
    options.Nx      (1,1) double {isinteger} = 128
    options.Ny      (1,1) double {isinteger} = 128
    options.GPU     (1,1) logical            = false
    options.Lx      (1,:) double             = 2*pi
    options.Ly      (1,:) double             = 2*pi
    options.dealias (1,1) double             = 0        % fraction of domain to dealias (e.g. 1/3)
    options.filter  (1,:) char               = 'none'   % 'none', 'exp'
    options.K_inner (1,1) double             = 2/3      % filter parameters
    options.K_outer (1,1) double             = 1        % filter parameters
    options.tol     (1,1) double             = 1e-15    % filter parameters
    options.order   (1,1) double             = 4        % filter parameters
    options.L       (:,:) double             = ones(128)
    options.N             function_handle    = @(u, grid, t) 0  % pass function using @(u, grid, t) function(u, grid, t)
    options.t       (1,1) double             = 0
    options.dt      (1,1) double             = 0.1
    options.iter    (1,1) double {isinteger} = 0
    options.stepper (1,:) char               = 'RK4'    % 'Euler', 'RK2', 'RK4'
    options.u       (:,:) double             = zeros(128)
end

grid = CreateGrid(Nx = options.Nx, Ny = options.Ny, GPU = options.GPU, Lx = options.Lx, ...
    Ly = options.Ly, dealias = options.dealias, filter = options.filter, K_inner = options.K_inner, ...
    K_outer = options.K_outer, tol = options.tol, order = options.order);

prob.grid = grid;

if options.GPU
    prob.u = gpuArray(options.u);
    prob.L = gpuArray(options.L);
else
    prob.u = options.u;
    prob.L = options.L;
end

prob.N = options.N;
prob.t = options.t;
prob.dt = options.dt;
prob.iter = options.iter;
prob.stepper = options.stepper;

end