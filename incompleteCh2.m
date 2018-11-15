function L2 = incompleteCh2(A)
% incomplete Cholesky factorisation with aij not eq 0
% N should be column length of A
L2 = tril(A);
N = length(A);
for i = 1:N
    for j = 1:N
        if A(i,j)~= 0
            for k = 1:N
                L2(k,k) = sqrt(L2(k,k));
                for i = k+1:N
                    L2(i,k) = L2(i,k)/L2(k,k);
                end
                for j = k+1:N
                    for i = j:N
                        L2(i,j) = L2(i,j) - L2(i,k)*L2(j,k);
                    end
                end
            end
        end
    end
end
end
