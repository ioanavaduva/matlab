function [V, W] = get_rk_basis(A, poles, V, W)
% Function that generates a Rational Krylov basis 
    
    nV = size(V);
    I = eye(size(A));
    
    W = (A + poles*I)\W; 
%     K = [zeros(length(W)), W];
    KK = orth(W);
    VW = [V, KK(:)]; 
  
%     s = svd(VW);
    
    V = orth(VW);
    if size(V) == nV
        warning('Space dimension does not increase');
    end
end

