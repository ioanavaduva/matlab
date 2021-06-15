function [shifts, it, pole_norm] = irka_shifts_dst(A, rhs, poles, tol)
    ch = 1; 
    it = 0;
    [Q,D]=eig(A); d = diag(D);
%     n = size(A, 1);
%     theta = ((n:-1:1)*pi/(n+1))';
%     d = 2+2*cos(theta);
    poles = sort(poles, 'descend');
    while ch > tol %&& it<=10 % while relative change in the shifts is >tol
        it = it+1;
        V = get_irka_basis_dst(Q,d, rhs, poles);
        Ap = V'*A*V; 
        shifts = sort(eig(Ap), 'descend'); %keyboard
        ch = norm(poles-shifts)/norm(poles); % c = max(abs((poles-shifts)./(eps+poles)));
        poles = shifts;
        pole_norm(it) = ch;
    end
end