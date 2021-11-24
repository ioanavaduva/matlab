function V = get_rk_basis_complex(A, poles, V)
% Function that generates a Rational Krylov basis 
    
    % Need to adjust pole set so that if there are complex conjugates one 
    % is removed from the pole set. 

    nV = size(V);
    I = eye(size(A));
    
    w = V(:, end); %keyboard
    w = (A + poles*I)\w; %keyboard
    
    wr = real(w); wi = imag(w);
    
    wr = wr - V*(V'*wr); %keyboard
    wr = wr - V*(V'*wr);%keyboard
    wr = wr/norm(wr); %keyboard
    V = [V, wr(:)]; 
    
    wi = wi - V*(V'*wi); %keyboard
    wi = wi - V*(V'*wi);%keyboard
    wi = wi/norm(wi); %keyboard
    V = [V, wi(:)];  

    if size(V) == nV
        warning('Space dimension does not increase');
    end
end

