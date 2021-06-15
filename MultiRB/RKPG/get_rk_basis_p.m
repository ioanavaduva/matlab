function V = get_rk_basis_p(A, Q, d, poles, V)
% Function that generates a Rational Krylov basis 
    
    nV = size(V);
%     I = eye(size(A));
   
    w = V(:, end); %keyboard
    w = Q*((Q'*w)./(d+poles));
    w = w - V*(V'*w); %keyboard
    w = w - V*(V'*w);%keyboard
    w = w/norm(w); %keyboard

    V = [V, w(:)]; 

    if size(V) == nV
        warning('Space dimension does not increase');
    end
end
