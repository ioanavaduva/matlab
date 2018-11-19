function L1 = incompleteCh3(A, TOL)
% incomplete Cholesky factorisation with |lij|<tol
% N should be column length of A
A = sparse(A);
L1 = tril(A);
N = length(A);
for k = 1:N
    L1(k,k) = sqrt(L1(k,k));
    L1(k+1:N,k) = L1(k+1:N,k)/L1(k,k);
    w = L1(k+1:N,k);
    w(abs(w)<TOL) = 0;
    L1(k+1:N,k) = w;
    
    for j = k+1:N
            L1(j:N,j) = L1(j:N,j) - L1(j:N,k)*L1(j,k);
    end
end

L1(abs(L1)<TOL) = 0;
% for i = 1:N 
%     for j = 1:N
%         if abs(L1(i,j)) < TOL
%             L1(i,j) = 0;
%         end
%     end
% end
% end

