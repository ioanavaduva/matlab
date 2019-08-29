function [x, maxit] = kronGenSylv(n, A, B, M1, M2, N1, N2, rhs, x0, tol, maxit)
    Nx = (- kron(N1', M1) - kron(N2', M2))*x0;
    for i = 1:maxit
        Mx = Nx + rhs;
        x = (kron(A, eye(n)) + kron(eye(n), B))\Mx;
        Nx = (- kron(N1', M1) - kron(N2', M2))*x;
        R = (kron(A, eye(n)) + kron(eye(n), B))*x - Nx - rhs;
        if norm(R) < tol
            maxit = i;
            return
        end
    end
end

