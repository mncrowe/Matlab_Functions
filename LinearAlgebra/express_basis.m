function c = express_basis(c, V)
% Expresses c in terms of the basis V = [e1 e2 ... en]
%
% - c: vector
% - V: matrix with basis vectors as columns

arguments
    c (:,1) double
    V (:,:) double
end

c = (V' * V) \ (V' * c);

end

