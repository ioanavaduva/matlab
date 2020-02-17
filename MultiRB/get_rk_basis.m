function [V, W] = get_rk_basis(A, poles, V, W)
% Function that generates a Rational Krylov basis 
    
    nV = size(V);
    I = eye(size(A));
    W
    W = (A + poles*I)\W 
    keyboard
    VW = [V, W(:)];  
  
    s = svd(VW);
    
    V = orth(VW);
    if size(V) == nV
        warning('Size of space does not increase');
    end
end

