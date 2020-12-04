function B_hat = tangential_dir(A, B, poles)
% truncated SVD of RHS
[U, ~, ~] = svd(B);
u = U(:, 1:length(poles));
% project A and B
A_r = u'*A*u;
B_r = u'*B;
% eigenvalue decomposition of A_r
[V, ~] = eig(A_r);
% tangential directions are columns of B_hat
B_hat = B_r'*V; %keyboard
end


