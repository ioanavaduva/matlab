function V = get_rk_basis(A, poles, V)
% Function that generates a Rational Krylov basis 
    
    nV = size(V);
    I = eye(size(A));
    
    w = V(:, end);
    w = (A + poles*I)\w; 
    w = w - V*(V'*w); w = w - V*(V'*w);
    w = w/norm(w);

    V = [V, w(:)]; 

    if size(V) == nV
        warning('Space dimension does not increase');
    end
end

