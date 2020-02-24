function V = get_rk_basis(A, poles, V)
% Function that generates a Rational Krylov basis 
    
    nV = size(V);
    I = eye(size(A));
    
    w = V(:, end); %keyboard
    w = (A + poles*I)\w; %keyboard
    w = w - V*(V'*w); %keyboard
    w = w - V*(V'*w);%keyboard
    w = w/norm(w); %keyboard

    V = [V, w(:)]; 

    if size(V) == nV
        warning('Space dimension does not increase');
    end
end

