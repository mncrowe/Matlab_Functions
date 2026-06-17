function [e, B, I] = orthog_space(v)
% Consider k vectors v_i in R^n. This function finds n-k vectors e_i which
% are orthogonal to the v_i and mutually orthogonal. The v_i are also
% orthogonalised in the process.
%
% Input:
% - v: k vectors, enter as v = [v_1, v_2, ..., v_k]
%
% Outputs:
% - e: basis vectors e = [e_1, ..., e_{N-k}]
% - B: full orthonormal basis of R^n
% - I: position of elements of e within B
%
% -------------------------------------------------------------------------
% Note: the v_i are extended to a full basis by including elements of the
% standard basis. The v_i only replace elements they have a component in.
% This basis is then converted to an orthonormal one using the Gram-Schmidt
% process.
% -------------------------------------------------------------------------

arguments
    v (:,:) double
end

s = size(v);
N = s(1);
k = s(2);
eps_err = 1e-6;

B = eye(N);
I = 1:N;

for i = 1:k
    j = 1;
    while length(I)>N-i
        if j>length(I); error('The v must be linerly independent.'); end
        if v(:,i)'*B(:,I(j)) > eps_err
            B(:,I(j)) = v(:,i);
            I(j) = [];
        end
        j = j+1;
    end
end

for j = 1:N
    for i = 1:j-1
        B(:,j) = B(:,j) - B(:,i)*(B(:,i)'*B(:,j))/(B(:,i)'*B(:,i));
    end
end

B = B./vecnorm(B,2,1);
e = B(:,I);

end