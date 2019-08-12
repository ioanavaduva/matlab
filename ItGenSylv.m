function [X, maxit] = ItGenSylv(A, B, N, M, C, X0, tol, maxit)
% N = [N1, N2, ...] and M = [M1, M2, ...]
    Nx = - cell2mat(M(1))*X0*cell2mat(N(1)) - cell2mat(M(2))*X0*cell2mat(N(2));
    Mx = Nx-C;
    for i = 1:maxit
        X = lyap(A, B, -Mx);
        R = A*X + X*B - Nx + C;
        if norm(R) < tol
            maxit = i;
            return
        end
        Nx = - cell2mat(M(1))*X*cell2mat(N(1)) - cell2mat(M(2))*X*cell2mat(N(2));
        Mx = Nx-C;
    end
end

