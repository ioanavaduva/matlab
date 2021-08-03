function V = get_irka_basis(A, rhs, poles)
    I = speye(size(A));
    V = zeros(length(A(:, 1)), length(poles));
    
    v0 = (A + poles(1)*I)\rhs; 
    V(:, 1) = v0/norm(v0);
    
    for i = 2:length(poles)
        w = V(:, i-1); 
        w = (A + poles(i)*I)\w; 
        w = w - V*(V'*w); 
        w = w - V*(V'*w);
        w = w/norm(w); 
        V(:, i) = w(:); 
    end
end