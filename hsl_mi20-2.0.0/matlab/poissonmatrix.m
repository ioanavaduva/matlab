function A = poissonmatrix(n)  
%n is the number of unknowns U
A = diag(4*ones(n^2, 1)) + diag (-1*ones(n^2-1, 1), 1) + diag (-1*ones(n^2-1, 1), -1) + diag (-1*ones(n^2-n, 1), n) + diag (-1*ones(n^2-n, 1), -n);
M = zeros(n^2, n^2);
for i = 1:n-1
    M(i*n, i*n+1) = 1;
    M(i*n+1, i*n) = 1;
end
A = A + M;
end

