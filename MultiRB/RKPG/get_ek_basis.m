function V = get_ek_basis(A, V)
% Function that generates a Rational Krylov basis 
    
    nV = size(V);
    if size(V, 2) < 2
        w = V(:,end); 
    else
        w = V(:, end-1);
    end
        
    if mod(size(V, 2), 2) == 0
        w = A\w; 
    else
        w = A * w; 
    end

    w = w - V*(V'*w); 
    w = w - V*(V'*w);
    w = w / norm(w);    
    
    V = [V, w(:)];
   
    if size(V) == nV
        warning('Space dimension does not increase');
    end
end
