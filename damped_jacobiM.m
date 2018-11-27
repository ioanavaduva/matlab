function x = damped_jacobiM(w, x0, A, b, TOL, maxit) 
M = sparse(diag(diag(A)));
N = sparse(-(tril(A, -1) + triu(A, 1)));
it = maxit;
for i = 1:it
    x = w*(M\(N*x0) + M\b) + (1-w)*x0;
    if norm(x-x0, inf) < TOL
        it = i;
        break;
    end
    x0 = x;
end
end

