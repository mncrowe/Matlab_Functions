function [Fp, xp, yp, coeffs, order] = poly_interp(F, p, x, y, options)
% Interpolates F to the point p where x (and y) correspond to the specified
% dimensions of F. This function can do 1D and 2D interpolation.
%
% Inputs:
% - F: muti-dimensional array with values to interpolate
% - p: point to interolate to, if using find_max = true, p is not used
% - x: coordinate of first interpolation dimension (default: 1:Nx)
% - y: corrdinate of second interpolation dimension, not used for 1D interpolation (default: 1:Ny)
% - options:
%   - n: number of data points in each direction used to fit coefficients (default: 3)
%   - m: polynomial order in each coordinate and maximum order (default: [2 2 2])
%   - find_max: true - find the maximum of F by finding the max of the interpolating polynomial (default: false)
%   - dims: dimension(s) to interpolate in (default: [1 2])
%   - closest_pts: true - use only the closest point to p when fitting coefficients (default: false)
%   - tol: tolerance for determining the maximum of F with 'fsolve' (default: 1e-8)
%   - debug: true - display information regarding identification of maximum (default: false)
%
% Outputs:
% - Fp: interpolated value of F at point p or at location of the maximum
% - xp: array containing the value of p(1) or the x position of the maximum
% - yp: array containing the value of p(2) or the y position of the maximum
% - coeffs: array containing the coefficients of the interpolating polynomial
% - order: the order of x (and y) in the term corresponding to each coefficient
%
% -------------------------------------------------------------------------
% Notes:
%
% This function can be used to interpolate a function F(x, y) to a given
% point p. It can also be used (by setting find_max = true) to determine
% the maximum value of F(x, y) by maximising the value of the interpolating
% polynomial around the region of the largest value of F. If find_max =
% true is used, the value of p can be anything.
%
% The parameter n determines the number of points in each direction which
% are used. This n x n grid is centred on the closest gridpoint to p. There
% must be a sufficient number of points to determine all % coefficients in
% the polynomial interpolation. If closest_pts = true, the minimum number
% of points is used such that the coefficients can be determined, these
% points are chosen to be those points which are closest to p and still lie
% with the n x n grid. If the unmber of points used to determine the
% coefficients exceeds the number of coefficients, linear regression is
% used to determine the best fit of the coefficients to the given points.
%
% The parameter m determines the degree of the interpolating polynomial in
% x and y. m is specified as an array with 2 or 3 entries. The first entry
% determines the degree in x, the second gives the degree in y, and the
% third gives the maximum order of each term. If only 2 entries are given,
% the third entry is set to the sum of the first 2. For example:
%   - m = [1 1 2] gives the polynomial: F = a + bx + cy + dxy
%   - m = [1 1 1] gives the polynomial: F = a + bx + cy
%   - m = [2 2 2] gives the polynomial: f = a + bx + cy + dx² + exy + fy²
%
% The dimensions of F to be interpolated through are specified by 'dims'.
% The default is dims = [1 2], this corresponds to 2D interpolation through
% the first two dimensions of F. The grid x should correspond to the
% dimension dims(1) and the grid y should correspond to the dimension 
% dims(2). This function may also be used to do 1D interpolation, in which
% case dims should be specified as a single integer and the value of y is
% not used (or required). Similarly, only the first value of p is used so p
% can (nd should) be specified as an integer.
%
% When using find_max = true, polynomials with m = [2 2 2] can be maximised
% analytically. For other values of m, a nonlinear root finding method is
% used. We use 'fsolve' which requires the MATLAB optimisation toolbox. The
% parameter 'tol' defines the tolerance for fsolve which is used to
% inentify when a root (corresponding to a maximum of F) has been found.
% The parameter 'debug' may be set to 'true' to display information
% regarding the identification of the maximum. In particular, the function
% will point out when the maximum occurs at a grid boundary and will
% display convergence information from the root finding step.
% -------------------------------------------------------------------------

arguments
    F                          double
    p                    (:,1) double
    x                    (:,1) double                         = (1:size(F,1))'
    y                    (1,:) double                         = 1:size(F,2)
    options.n            (1,:) double {mustBeInteger}         = 3
    options.m            (1,:) double {mustBeInteger}         = [2 2 2]
    options.find_max     (1,1) logical                        = false
    options.dims         (1,:) double {mustBeInteger,MaxDims} = [1 2]
    options.closest_pts  (1,1) logical                        = false
    options.tol          (1,1) double                         = 1e-8
    options.debug        (1,1) logical                        = false
end

% Set options for fsolve is using find_max = true:

if options.find_max
    if options.debug
        fsolve_opts = optimoptions('fsolve', 'TolFun', ...
            options.tol, 'SpecifyObjectiveGradient', true);
    else
        fsolve_opts = optimoptions('fsolve', 'Display', 'none', 'TolFun', ...
            options.tol, 'SpecifyObjectiveGradient', true);
    end
end

% Require polynomial to be quadratic or higher if using find_max = true:

if options.find_max && min(options.m) < 2
    eid = 'FindMax:InterpError';
    msg = 'Using find_max = true requires the interpolating polynomial to be quadratic or higher.';
    error(eid, msg)
end

% Set to a 1D problem if F has only 1 dimension:

if length(F) == numel(F)
    options.dims = 1;
end

% Rearrange interpolation dimensions to the start of F:

I = 1:length(size(F));
I(options.dims) = [];

F = permute(F, [options.dims I]);

% Move all remaining dimensions to a single dimension:

s = size(F);

if length(options.dims) == 1
    F = reshape(F, s(1), []);
else
    F = reshape(F, s(1), s(2), []);
end

% Deal with 1D and 2D cases seperately:

if length(options.dims) == 1 % 1D case:

    % Set 1D m and n using first entry of m and n only:

    options.m = options.m(1);
    options.n = options.n(1);

    % Check length of x in consistent with F:

    if length(x) ~= size(F, 1)
        eid = 'Grid:InconsistentLength';
        msg = ['The length of x is not consistent with dimension ' num2str(options.dims) ' of F.'];
        error(eid, msg)
    end

    % Ensure sufficient points are used to calculate all polynomial coefficients:

    if  options.n < options.m + 1 
        eid = 'Coefficients:InsufficientData';
        msg = ['n = ' num2str(options.n) ' data points is insufficient to determine the ' ...
            num2str(options.m + 1) ' coefficients for a polynomial of degree m = ' num2str(options.m) '.'];
        error(eid, msg)
    end

    if options.find_max % Find maximum using polynomial interpolation (1D):

        % Define arrays:

        coeffs = zeros(options.m + 1, size(F, 2));
        Fp = zeros(1, size(F, 2));
        xp = Fp;

        % Loop through second dimension of F, finding maximum over dim 1:
        
        for i = 1:size(F, 2)

            % Identify location of maximum of F:
    
            [~, ix] = max(F(:, i));

            xp(i) = x(ix);

            if ismember(ix, [1 length(x)]) % Maximum occurs at an edge:

                % Warn if the maximum occurs at the edge of the domain:

                if options.debug
                    disp(['Maximum occurs at edge (x = ' num2str(x(ix)) ').'])
                end

                % Set maximum to value of F at domain edge in this case:

                Fp(1, i) = F(ix, i);

            else % Maximum occurs in the interior:

                % Get n points closest to maximum and fit coefficients:
    
                [~, ix] = sort(abs(x - xp(i)));
                ix = ix(1:options.n);
    
                coeffs(:, i) = GetCoeff_1D(x(ix), F(ix, i), options.m);

                % Set Fp and xp using fsolve to find maximum of interpolating polynomial:

                [Fp(i), xp(i)] = FindMax_1D(coeffs(:, i), xp(i), fsolve_opts, options.debug);

            end

        end

    else % Find value of F at p using polynomial interpolation (1D):

        % Get n closest points to p:

        [~, ix] = sort(abs(x - p(1)));
        ix = ix(1:options.n);

        % Define empty arrays:

        coeffs = zeros(options.m + 1, size(F, 2));
        Fp = zeros(1, size(F, 2));
        xp = p * ones(1, size(F, 2));

        % Loop through second dimension of F, finding coefficients and interpolating to F(p):
        
        for i = 1:size(F, 2)

            coeffs(:, i) = GetCoeff_1D(x(ix), F(ix, i), options.m);
  
            Fp(1, i) = sum(coeffs(:, i) .* (p .^ (0:options.m))');

        end

    end

    % Reshape output arrays to match initial number of dimensions:

    Fp = reshape(Fp, [1  s(2:end)]);
    xp = reshape(xp, [1  s(2:end)]);
    yp = 0*xp;

    % Output the order of x correspoding to each coefficient:

    order = 0:options.m;

else % 2D case:

    % Check length of x and y are consistent with size of F:

    if length(x) ~= size(F, 1)
        eid = 'Grid:InconsistentLength';
        msg = ['The length of x is not consistent with dimension ' num2str(options.dims(1)) ' of F.'];
        error(eid, msg)
    end

    if length(y) ~= size(F, 2)
        eid = 'Grid:InconsistentLength';
        msg = ['The length of y is not consistent with dimension ' num2str(options.dims(2)) ' of F.'];
        error(eid, msg)
    end

    % If only 1 m value is given, set order in both x and y directions to m:

    if length(options.m) == 1
        options.m = [options.m options.m];
    end

    % If no maximum term order is given, take the maximum as m(1) + m(2):

    if length(options.m) == 2
        options.m = [options.m options.m(1)+options.m(2)];
    end

    % of only one n value s given, take n points in both x and y directions:

    if length(options.n) == 1
        options.n = [options.n options.n];
    end

    % Determine number of coefficients required using options.m input:

    N_coeff = (options.m(1) + 1) * (options.m(2) + 1) - (options.m(1) + options.m(2) - ...
        min(options.m(1) + options.m(2), options.m(3))) * (options.m(1) + options.m(2) ...
        - min(options.m(1) + options.m(2), options.m(3)) + 1) / 2;

    % Ensure sufficient points are used to calculate all polynomial coefficients:

    if options.n(1) * options.n(2) < N_coeff
        eid = 'Coefficients:InsufficientData';
        msg = ['n(1)*n(2) = ' num2str(options.n(1) * options.n(2)) ' data points is insufficient to determine the ' ...
            num2str(N_coeff) ' coefficients for a polynomial of degree (m(1), m(2)) = (' num2str(options.m(1)) ', ' ...
            num2str(options.m(2)) ') with a maximum order per term of m(3) = ' num2str(options.m(3)) '.'];
        error(eid, msg)
    end

    if options.find_max % Find maximum using polynomial interpolation (2D):

        % Define empty arrays:

        Fp = zeros(1, 1, size(F, 3));
        xp = Fp;
        yp = Fp;
        coeffs = zeros(N_coeff, size(F, 3));
        order = zeros(N_coeff, 2);

        % Loop through third dimension of F, finding maximum over dims 1 & 2:

        for i = 1:size(F, 3)

            % Identify location of maximum of F:
    
            [Fm, ix] = max(F(:, :, i));
            [~, iy] = max(Fm);
            ix = ix(iy);

            xp(i) = x(ix);
            yp(i) = y(iy);

            if ismember(ix, [1 length(x)])

                if ismember(iy, [1 length(y)]) % Maximum occurs at a corner:

                    if options.debug
                        disp(['Maximum occurs at corner (x, y) = (' num2str(x(ix)) ', ' num2str(y(iy)) ').'])
                    end

                    % Set maximum to value of F at domain edge in this case:

                    Fp(1, 1, i) = F(ix, iy, i);

                else % Maximum occurs along an x boundary:

                    if options.debug
                        disp(['Maximum occurs at edge (x = ' num2str(x(ix)) ').'])
                    end

                    % Get n(2) points closest to maximum and fit coefficients:

                    [~, iy] = sort(abs(y - yp(i)));
                    iy = iy(1:options.n(2));

                    coeffs(1:options.m(2)+1, i) = GetCoeff_1D(y(iy)', F(ix, iy, i)', options.m(2));
                    order(1:options.m(2)+1, 2) = 0:options.m(2);

                    % Set Fp and yp using fsolve to find maximum of interpolating polynomial:

                    [Fp(i), yp(i)] = FindMax_1D(coeffs(1:options.m(2)+1, i), yp(i), fsolve_opts, options.debug);

                end

            else

                if ismember(iy, [1 length(y)]) % Maximum occurs along a y boundary:

                    if options.debug
                        disp(['Maximum occurs at edge (y = ' num2str(y(iy)) ').'])
                    end
                    
                    % Get n(1) points closest to maximum and fit coefficients:

                    [ ~, ix] = sort(abs(x - xp(i)));
                    ix = ix(1:options.n(1));

                    coeffs(1:options.m(1)+1, i) = GetCoeff_1D(x(ix), F(ix, iy, i), options.m(1));
                    order(1:options.m(1)+1, 2) = 0:options.m(1);

                    % Set Fp and yp using fsolve to find maximum of interpolating polynomial:

                    [Fp(i), xp(i)] = FindMax_1D(coeffs(1:options.m(1)+1, i), xp(i), fsolve_opts, options.debug);

                else % Maximum occurs in the interior:

                    % Get closest x and y points to xp(i) and yp(i):

                    [~, ix] = sort(abs(x - xp(i)));
                    ix = ix(1:options.n(1));
            
                    [~, iy] = sort(abs(y - yp(i)));
                    iy = iy(1:options.n(2));
            
                    % Reshape x and y grids around p into vectors:
            
                    xm = reshape(x(ix) * ones(1, options.n(2)), [], 1);
                    ym = reshape(ones(options.n(1), 1) * y(iy), [], 1);

                    % Find coefficients:

                    fm = reshape(F(ix, iy, i), [], 1);
    
                    [coeffs(:, i), order] = GetCoeff_2D(xm, ym, fm, options.m, options.closest_pts, N_coeff, p);

                    % Set xp, yp and Fp using fsolve to find maximum of F:

                    [Fp(i), xp(i), yp(i)] = FindMax_2D(coeffs(:, i), order, xp(i), yp(i), fsolve_opts, options.debug, options.m);

                end

            end

        end

    else % Find value of F at p using polynomial interpolation (2D):

        % Get closest x and y points to p:

        [~, ix] = sort(abs(x - p(1)));
        ix = ix(1:options.n(1));

        [~, iy] = sort(abs(y - p(2)));
        iy = iy(1:options.n(2));

        % Reshape x and y grids around p into vectors:

        xm = reshape(x(ix) * ones(1, options.n(2)), [], 1);
        ym = reshape(ones(options.n(1), 1) * y(iy), [], 1);

        % Set empty arrays:

        coeffs = zeros(N_coeff, size(F, 3));
        Fp = zeros(1, 1, size(F, 3));

        % Set xp and yp to p values:

        xp = p(1) * ones(1, 1, size(F, 3));
        yp = p(2) * ones(1, 1, size(F, 3));

        % Loop through third dimension of F, finding coefficients and interpolating to F(p):
        
        for i = 1:size(F, 3)

            fm = reshape(F(ix, iy, i), [], 1);
    
            [coeffs(:, i), order] = GetCoeff_2D(xm, ym, fm, options.m, options.closest_pts, N_coeff, p);

            Fp(1, 1, i) = sum(coeffs(:, i) .* p(1) .^ order(:, 1) .* p(2) .^ order(:, 2));

        end

    end

    % Reshape output arrays to match initial number of dimensions:

    Fp = reshape(Fp, [1 1 s(3:end)]);
    xp = reshape(xp, [1 1 s(3:end)]);
    yp = reshape(yp, [1 1 s(3:end)]);

end

% Permute dimensions to match order of initial input array:

Fp = ipermute(Fp, [options.dims I]);
xp = ipermute(xp, [options.dims I]);
yp = ipermute(yp, [options.dims I]);

end

function MaxDims(dims)
    if length(dims) > 2
        eid = 'Dims:ExceedMax';
        msg = 'Polynomial interpolation can be done in at most 2 dimensions.';
        error(eid, msg)
    end
end

function C = GetCoeff_1D(x, f, m)
    M = x.^(0:m);
    C = M \ f;
end

function [C, Order] = GetCoeff_2D(x, y, f, m, closest_pts, N_coeff, p)

    Order = zeros(1 + m(1), 1 + m(2), 3); % Matrix containing exponents for x and y

    Order(:, :, 1) = (0:m(1))' * ones(1, 1 + m(2));
    Order(:, :, 2) = ones(1 + m(1), 1) * (0:m(2));
    Order(:, :, 3) = Order(:, :, 1) + Order(:, :, 2);

    Order = reshape(Order, [], 3);

    if closest_pts  % if cloests_pts = true, take only the points N_coeffs points closest to p
        [~, i] = sort((x - p(1)).^2 + (y - p(2)).^2);
        x = x(i(1:N_coeff));
        y = y(i(1:N_coeff));
        f = f(i(1:N_coeff));
    end
    
    M = x .^ (Order(:, 1)') .* y .^ (Order(:, 2)'); % Matrix containing fit data
    
    M(:, Order(:,3) > m(3)) = [];   % Remove columns with total order greater than m(3)
    Order(Order(:, 3) > m(3), :) = [];
    Order = Order(:, 1:2);
    
    C = M \ f;  % solve for coefficients using inversion or linear regression

end

function [F, J] = FJ_1D(x, C, m)
    i1 = 0:m-1;
    i2 = 0:m-2;

    C1 = (1:m) .* C(2:end)';
    C2 = (2:m) .* (1:m-1) .* C(3:end)';

    F = sum(C1 .* (x .^ i1));
    J = sum(C2 .* (x .^ i2));
end

function [F, J] = FJ_2D(x, y, C, Order)
    C1 = (Order(:, 1) .* C)';
    C2 = (Order(:, 2) .* C)';

    O11 = max(Order(:, 1)'-1, 0);
    O12 = Order(:, 2)';
    O21 = Order(:, 1)';
    O22 = max(Order(:, 2)'-1, 0);

    O111 = max(Order(:, 1)'-2, 0);
    O222 = max(Order(:, 2)'-2, 0);

    Fx = sum(C1 .* x .^ O11 .* y .^ O12);
    Fy = sum(C2 .* x .^ O21 .* y .^ O22);

    C11 = C1 .* O11;
    C12 = C1 .* O12; % = C2 .* O21
    C22 = C2 .* O22;

    Fxx = sum(C11 .* x .^ O111 .* y .^ O12);
    Fxy = sum(C12 .* x .^ O11 .* y .^ O22);
    Fyy = sum(C22 .* x .^ O21 .* y .^ O222);

    F = [Fx; Fy];
    J = [Fxx Fxy; Fxy Fyy];
end

function [Fp, xp] = FindMax_1D(coeffs, xp, fsolve_opts, debug)

    m = length(coeffs) - 1;             % Get polynomial order

    if m == 2   % Deal with quadratics analytically

        if debug; disp('Finding maximum analytically.'); end

        xp = -coeffs(2) / (2 * coeffs(3));
        Fp = coeffs(1) - coeffs(2)^2 / (2 * coeffs(3));

    else

        if debug; disp('Finding maximum using fsolve:'); end

        fJ = @(x) FJ_1D(x, coeffs, m);      % Define df/dx and d^2/dx^2 for fsolve

        xp = fsolve(fJ, xp, fsolve_opts);   % Find position of maximum
        Fp = sum(coeffs .* (xp .^ (0:m))'); % Evaluate interpolating polynomial at xp

    end
end

function [Fp, xp, yp] = FindMax_2D(coeffs, order, xp, yp, fsolve_opts, debug, m)

    if m(1) == 2 && m(2) == 2 && m(3) == 2 % Deal with quadratics analytically

        if debug; disp('Finding maximum analytically.'); end
        
        xp = (-2*coeffs(6)*coeffs(2) + coeffs(5)*coeffs(4)) / (4*coeffs(3)*coeffs(6) - coeffs(5)^2);
        yp = (coeffs(5)*coeffs(2) - 2*coeffs(3)*coeffs(4)) / (4*coeffs(3)*coeffs(6) - coeffs(5)^2);
        Fp = coeffs(1) + coeffs(2)*xp + coeffs(3)*xp^2 + coeffs(4)*yp + coeffs(5)*xp*yp + coeffs(6)*yp^2;

    else

        if debug; disp('Finding maximum using fsolve:'); end

        fJ = @(x) FJ_2D(x(1), x(2), coeffs, order); % Define Grad(F) and J[f] for fsolve

        p = fsolve(fJ, [xp; yp], fsolve_opts);                      % find position of maximum
        xp = p(1);                                                  % find xp
        yp = p(2);                                                  % find yp
        Fp = sum(coeffs .* xp .^ order(:, 1) .* yp .^ order(:, 2)); % Evaluate interpolating polynomial here

    end

end