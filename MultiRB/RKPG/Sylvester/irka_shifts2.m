function [shifts1, shifts2, it] = irka_shifts2(A, B, rhs1, rhs2, poles1, poles2, tol)
    ch1 = 1; 
    ch2 = 1;
    
    it = 0;
    A = sparse(A);
    B = sparse(B);
    
    poles1 = sort(poles1, 'descend');
    poles2 = sort(poles2, 'descend');
    
    while ch1 > tol && ch2 > tol %&& it<=10 % while relative change in the shifts is >tol
        it = it+1;
        V = get_irka_basis(A, rhs1, poles1);
        W = get_irka_basis(B, rhs2, poles2);
        
        Ap = V'*A*V; 
        Bp = W'*B*W;
        
        shifts1 = sort(eig(Ap), 'descend'); 
        shifts2 = sort(eig(Bp), 'descend');
        
        ch1 = norm(poles1-shifts2)/norm(poles1); % c = max(abs((poles-shifts)./(eps+poles)));
        ch2 = norm(poles2-shifts1)/norm(poles2); 
        
        poles1 = shifts2;
        poles2 = shifts1;

    end
end