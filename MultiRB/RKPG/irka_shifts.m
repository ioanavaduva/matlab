function [shifts, it] = irka_shifts(A, rhs, poles, tol)
    ch = 1; 
    it = 0;
    while ch >tol
        it = it+1;
        I = eye(size(A));
        v0 = (A + poles(1)*I)\rhs;
        V = v0;
        for i = 2:length(poles)
            V = rk_basis_irka(A, poles(i), V);
        end

        Ap = V'*A*V; 
        rhsp = V'*rhs;
        [Q, D] = eig(Ap);
        shifts = diag(D);
        new_rhs = rhsp'*Q;
        ch = norm(poles-shifts)/norm(poles);
        poles = shifts;
    end
end