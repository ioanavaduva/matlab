function [shifts, it] = irka_shifts_block(A, B, B_hat, poles, tol)
    ch = 1; 
    it = 0;
    poles = sort(poles, 'descend');
    while ch > tol % while relative change in the shifts is >tol
        it = it+1; %keyboard
        V = get_irka_block_basis(A, B, B_hat, poles);
        Ap = V'*A*V; 
        [Q, s] = eig(Ap);
        shifts = sort(diag(s), 'descend'); %keyboard
        ch = norm(poles-shifts)/norm(poles); %keyboard % c = max(abs((poles-shifts)./(eps+poles)));
        poles = shifts; 
        Bp = V'*B;
        B_hat = Bp'*Q;
    end
end