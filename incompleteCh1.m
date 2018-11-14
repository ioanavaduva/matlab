function L = incompleteCh1(A, TOL)
% incomplete Cholesky factorisation with |lij|<tol
% N should be column length of A
% Checked that Ch*Ch' = A
L = tril(A);
N = length(A);
for k = 1:N
    L(k,k) = sqrt(L(k,k));
    for i = k+1:N
        L(i,k) = L(i,k)/L(k,k);
    end
    for j = k+1:N
        for i = j:N
            L(i,j) = L(i,j) - L(i,k)*L(j,k);
        end
    end
end
for i = 1:N 
    for j = 1:N
        if abs(L(i,j)) < TOL
            L(i,j) = 0;
        end
    end
end

