function [X, maxit] = ItGenSylv(A, B, N, M, C, X0, tol, maxit)
% N = [N1, N2, ...] and M = [M1, M2, ...]
    Nx = - cell2mat(M(1))*X0*cell2mat(N(1)) - cell2mat(M(2))*X0*cell2mat(N(2));
    for i = 1:maxit
        Mx = Nx-C;
        X = lyap(A, B, -Mx);
        Nx = - cell2mat(M(1))*X*cell2mat(N(1)) - cell2mat(M(2))*X*cell2mat(N(2));
        R = A*X + X*B - Nx + C;
        norm(R)
        if norm(R) < tol
            maxit = i;
            return
        end
    end
end

