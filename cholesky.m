function Ch = cholesky(A)
% Cholesky factorisation
% N should be column length of A
% Must work with lower part of A to obtain a lower triangular matrix Ch -
% this is the same as built in chol(A, 'lower')
% Checked that Ch*Ch' = A
Ch = tril(A);
N = length(A);
for k = 1:N
    Ch(k,k) = sqrt(Ch(k,k));
    for i = k+1:N
        Ch(i,k) = Ch(i,k)/Ch(k,k);
    end
    for j = k+1:N
        for i = j:N
            Ch(i,j) = Ch(i,j) - Ch(i,k)*Ch(j,k);
        end
    end
end

