function L1 = incompleteCh1(A, TOL)
% incomplete Cholesky factorisation with |lij|<tol
% N should be column length of A
L1 = tril(A);
N = length(A);
for k = 1:N
    L1(k,k) = sqrt(L1(k,k));
    for i = k+1:N
        L1(i,k) = L1(i,k)/L1(k,k);
    end
    for j = k+1:N
        for i = j:N
            L1(i,j) = L1(i,j) - L1(i,k)*L1(j,k);
        end
    end
end
for i = 1:N 
    for j = 1:N
        if abs(L1(i,j)) < TOL
            L1(i,j) = 0;
        end
    end
end
end

