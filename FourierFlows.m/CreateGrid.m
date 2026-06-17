function grid = CreateGrid(options)

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
end

if numel(options.Lx) == 1
    options.Lx = options.Lx * [-0.5 0.5];
end

if numel(options.Ly) == 1
    options.Ly = options.Ly * [-0.5 0.5];
end

dx = (options.Lx(2) - options.Lx(1)) / options.Nx;
dy = (options.Ly(2) - options.Ly(1)) / options.Ny;

x = (options.Lx(1):dx:(options.Lx(2)-dx))';
y = options.Ly(1):dy:(options.Ly(2)-dy);

k = 2 * pi / (options.Lx(2) - options.Lx(1)) * [0:ceil(options.Nx/2-1) -floor(options.Nx/2):-1]';
l = 2 * pi / (options.Ly(2) - options.Ly(1)) * [0:ceil(options.Ny/2-1) -floor(options.Ny/2):-1];

%k = pi / dx * (-1:2/options.Nx:(1-2/options.Nx))';
%l = pi / dy * (-1:2/options.Ny:(1-2/options.Ny));

switch options.filter
    case 'none'
        filter = 1.0 * (abs(k) <= (1-options.dealias) * pi / dx) .* (abs(l) <= (1-options.dealias) * pi / dy);
    case 'exp'
        K = sqrt((k * dx / pi).^2 + (l * dy / pi).^2);
        decay = -log(options.tol) / (options.K_outer - options.K_inner) ^ options.order;
        filter = exp(- decay * (K - options.K_inner) .^ options.order);
        filter(K < options.K_inner) = 1;
    otherwise
        error('Filter:UnrecognisedValue', 'Filter type is not recognised.')
end

grid.GPU = options.GPU;

grid.x = x;
grid.y = y;

K2_inv = 1 ./ (k.^2 + l.^2);
K2_inv(isinf(K2_inv)) = 0;

if grid.GPU
    grid.k = gpuArray(k);
    grid.l = gpuArray(l);
    grid.K2 = gpuArray(k.^2 + l.^2);
    grid.K2_inv = gpuArray(K2_inv);
    grid.filter = gpuArray(filter);
else
    grid.k = k;
    grid.l = l;
    grid.K2 = k.^2 + l.^2;
    grid.K2_inv = K2_inv;
    grid.filter = filter;
end

end