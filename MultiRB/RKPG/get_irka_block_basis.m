function V = get_irka_block_basis(A, rhs, poles)
    V = zeros(length(A(:, 1)), length(poles)*length(rhs(1, :)));
    I = speye(size(A));
    k = 1; 
    for i = 1:length(poles)
        for j = 1:length(rhs(1, :))
            v = (A + poles(i)*I)\rhs(:, j);
            v = v - V*(V'*v); 
            v = v - V*(V'*v);
            v = v/norm(v);
            V(:, k) = v(:);
            k = k+1;
        end
    end
end