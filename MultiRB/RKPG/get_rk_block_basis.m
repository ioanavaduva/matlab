function V = get_rk_block_basis(A, poles, V, col)
% Function that generates a Rational Krylov basis 
    
    nV = size(V);
    I = eye(size(A));

    w = (A + poles*I)\V(:,col);
    w = w - V*(V'*w); %keyboard
    w = w - V*(V'*w);%keyboard
    w = w/norm(w); %keyboard

    V = [V, w(:)]; %keyboard

    if size(V) == nV
        warning('Space dimension does not increase');
    end
end