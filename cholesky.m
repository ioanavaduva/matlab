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

% [~, p] = chol(matrix) gives p=0 is matrix is positive definite and any
% other value for p if it isn't
% my cholesky function gives correct matrix A when doing Ch*Ch' for
% positive definite and semidefinite matrices, otherwise gives irrational
% entries in the Ch matrix and Ch*Ch' doesn't give back original matrix