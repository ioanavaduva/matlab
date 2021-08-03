function V = get_irka_basis(A, rhs, poles)
% the basis will be given by span{(A+s1I)^{-1}b, ...,(A+skI)^{-1}b}
    I = speye(size(A));
    V = zeros(length(A(:, 1)), length(poles));
    
    v0 = (A + poles(1)*I)\rhs; %compute first column to add to basis
    V(:, 1) = v0/norm(v0);
    
    for i = 2:length(poles)
        w = V(:, i-1); %keyboard
        w = (A + poles(i)*I)\w; 
        w = w - V*(V'*w); %keyboard
        w = w - V*(V'*w);%keyboard
        w = w/norm(w); %keyboard
        V(:, i) = w(:); 
    end
end