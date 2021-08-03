function [shifts, it] = irka_shifts(A, rhs, poles, tol)
    ch = 1; 
    it = 0;
%     A = sparse(A);
    poles = sort(poles, 'descend');
    
    while ch > tol 
        it = it+1;
        V = get_irka_basis(A, rhs, poles);
        Ap = V'*A*V; 
        shifts = sort(eig(Ap), 'descend'); 
        ch = norm(poles-shifts)/norm(poles); 
        poles = shifts;

    end
end