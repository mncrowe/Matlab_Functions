function [phi,lambda] = POD(F, n)
% Performs a proper orthogonal decomposition of a field F
% 
% Inputs:
% - F: input field, size = [Nx,...,Nt]
% - n: number of POD modes to output
%
% Outputs:
% - phi: n POV modes, size = [Nx,...,N]
% - lambda: variance associated with phi, size = [N,1]

arguments
    F       double
    n (1,1) double {mustBeInteger}
end

s = size(F);

F = reshape(F, [prod(s(1:end-1)) s(end)]);
C = 1 / (s(end)-1) * F * (F');    % covariance matrix

[phi, lambda] = eigs(C, n);

lambda = sum(lambda)';
phi = reshape(phi, [s(1:end-1) numel(phi) / prod(s(1:end-1))]);

end

