% poisson matrix using kronecker tensor product
function AA = kr_pois(n)
T = diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1);
I = speye(n);
T = sparse(T);
AA = kron(T, I) + kron(I, T);
end

% resulting AA is sparse