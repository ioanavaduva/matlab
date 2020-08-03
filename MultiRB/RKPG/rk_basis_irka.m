function V = rk_basis_irka(A, poles, V)
    I = eye(size(A));
    
    w = V(:, end); %keyboard
    w = (A + poles*I)\w; %keyboard
    w = w - V*(V'*w); %keyboard
    w = w - V*(V'*w);%keyboard
    w = w/norm(w); %keyboard
    V = [V, w(:)]; 
end