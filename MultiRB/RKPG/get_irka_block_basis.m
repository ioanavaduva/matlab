function V = get_irka_block_basis(A, B, B_hat, poles)
I = speye(size(A));
V = (A + poles(1)*I)\(B*B_hat(:, 1)); 
V = V/norm(V); %[Q, ~] = qr(V); %should give orth basis for V
for i = 2:length(poles)
    v = (A + poles(i)*I)\(B*B_hat(:, i));
    v = v - V*(V'*v); %use Q instead of V        
    v = v - V*(V'*v);
    v = v/norm(v);
    V(:, i) = v(:);
end
end
    