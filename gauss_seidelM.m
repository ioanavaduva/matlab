function [it, x] = gauss_seidelM(x0, A, b, TOL, maxit)
%x0 is the initial guess, A is the matrix corresponding to the system, 
%b is a known vector, TOL = tolerance and maxit is the maximum number
%of iterations
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);
M = D+L;
N = -U;
it = maxit;
for i = 1:it
    x = M\(N*x0) + M\b;
    if norm(x-x0, inf) < TOL
        it = i;
        break;
    end
    x0 = x;
end
end


