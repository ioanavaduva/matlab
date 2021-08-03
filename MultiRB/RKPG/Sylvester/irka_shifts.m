function [shifts, it] = irka_shifts(A, rhs, poles, tol)
    ch = 1; 
    it = 0;
    A = sparse(A);
    poles = sort(poles, 'descend');
    
    while ch > tol %&& it<=10 % while relative change in the shifts is >tol
        it = it+1;
        V = get_irka_basis(A, rhs, poles);
        Ap = V'*A*V; 
        shifts = sort(eig(Ap), 'descend'); %keyboard
        ch = norm(poles-shifts)/norm(poles); % c = max(abs((poles-shifts)./(eps+poles)));
        poles = shifts;

    end
end