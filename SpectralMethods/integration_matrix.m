function I = integration_matrix(M, i, C)
% Defines the integration matrix which inverts the system My = f subject to 
% boundary condition Cy = d for row vector C and scalar d. This gives 
% solution y = If' where f' is f with the i^th row replaced by d.
%
% - M: differentiation matrix (n x n)
% - i: index of row to replace with boundary condition (default: 1)
% - C: boundary condition row (1 x n) (default: [0, .., 1,.., 0], 1 in position i)
%
% -------------------------------------------------------------------------
% Note: Since differentiation maps constant vectors to 0, the matrix M must
%       have a zero eigenvalue (with eigenvector 1) and hence a zero
%       determinant. This represents the need to impose an integration
%       constant.
% -------------------------------------------------------------------------

arguments
    M  (:,:) double
    i  (1,1) double {isinteger} = 1
    C  (1,:) double             = (1:length(M))==i
end

M(i, :) = C;
I = M^-1;

end

