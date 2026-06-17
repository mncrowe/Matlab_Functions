function [a, G, phi] = POD_Coeffs(F, phi)
% Calculates the POV coefficients for a field F and POV modes phi
%
% Inputs:
% - F: input field, should be size = [Nx,Ny,...,Nt] for Nt timesteps
% - phi: array or scalar;
% 	- array: POD modes, should be size = [Nx,Ny,...,N] for N modes
% 	- scalar: N, number of POV modes (default: 10)
%
% Outputs:
% - a: coefficient array, size = [N,Nt]
% - G: reconstructed F using given POD modes, size = [Nx,Ny,...,Nt]
% - phi: N POD modes, size = [Nx,Ny,...,N]
%
% Notes: Modes are not required to be normalised. Field F can be
% reconstructed by F = Sum_i [phi_i a_i] with an error corresponding to
% the neglected modes in phi.

arguments
    F   double
    phi double = 10
end

% If POD modes are not given, calculate them:

if numel(phi) == 1; phi = POD(F, phi); end

% Define output:

sp = size(phi);
sf = size(F);
a = zeros(sp(end), sf(end));

% Reshape F and phi into matrices:

F = reshape(F, [prod(sf(1:end-1)) sf(end)]);
phi = reshape(phi, [prod(sp(1:end-1)) sp(end)]);

% Project F onto POD modes:

for in = 1:sp(end)
    a(in, :) = (F') * phi(:,in) / ((phi(:,in)') * phi(:,in));
end

% Define G and reshape output back to original dimensions:

G = reshape(phi*a, sf);
phi = reshape(phi, sp);

end

