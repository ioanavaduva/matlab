function V = get_poly_basis(A, V)
% Function that generates a Rational Krylov basis 
    
    nV = size(V);
    
    w = V(:,end);
    w = A * w;                        
    w = w - V*(V'*w); 
    w = w - V*(V'*w);
    w = w / norm(w);    
    
    V = [V, w(:)];
   
    if size(V) == nV
        warning('Space dimension does not increase');
    end
end
