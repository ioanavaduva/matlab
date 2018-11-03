function [iter, x] = gauss_seidel(x0, A, b, TOL, N)
% x0 is initial guess, TOL is tolerance, N is max no of iterations
n = length(b);
k = 1;
iter = N;
x = x0;
xold = x0;
for k = 1:N
    for i = 1:n
        x(i) = (b(i) - sum(A(i, 1:i-1).*x(1:i-1)') - sum(A(i, i+1:n).*x(i+1:n)'))/A(i,i);
    end
    if norm(x-xold,inf) < TOL
        iter = k;
        break;
    end
    xold = x;
end
end
