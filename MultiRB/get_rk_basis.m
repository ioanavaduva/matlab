function [V, W] = get_rk_basis(A, poles, V, W)
% Function that generates a Rational Krylov basis 
    
    nV = size(V);
    I = eye(size(A));
    W = (A + poles*I)\W;
    VW = [V, W(:)];
    if size(VW) == nV
        warning('Size of space does not increase');
    end
    V = orth(VW);
    
end

