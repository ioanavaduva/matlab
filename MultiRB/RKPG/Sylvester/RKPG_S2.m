 function [X1, X2, final_err, vec_res, it] = RKPG_S2(M, N, rhs1, rhs2, poles1, poles2, tol, maxit)

    i = 0; 
    res = 1; 
    it = 0; % the iterations also act as space dimension
    
    v0 = rhs1/norm(rhs1);
    V = v0; 
    
    w0 = rhs2/norm(rhs2);
    W = w0;
    
    vec_res(1) = res;

    fprintf(' no.its  residual \n')
    
    while (res > tol && it < maxit)
        it = it + 1;
        i = i+1; 
        if i > length(poles1) % cycle through poles          
            i = 1; 
        end
        
        % choose basis 
        V = get_rk_basis(M, poles1(i), V); 
        W = get_rk_basis(N', poles2(i), W); 


        % project matrix A and rhs1/2
        Mp = V'*M*V; 
        Np = W'*N*W;
        rhs1p = V'*rhs1;
        rhs2p = W'*rhs2;%keyboard;
        
        % solve projected problem
         Y = lyap(-Mp, -Np, rhs1p*rhs2p'); 

        % obtain low-rank factors X_1 & X_2 using SVD
        [uu, ss, vv] = svd(Y, 0);
        X1 = V*uu*ss;
        X2 = vv'*W';
        
%         X = V*Y*W';
        
        % compute residual
        res = norm(M*X1*X2 + X1*X2*N - rhs1*rhs2', 'fro')/norm(rhs1*rhs2', 'fro');
%         res = norm(M*X + X*N - rhs1*rhs2', 'fro')/norm(rhs1*rhs2', 'fro');
        
        % Print details to screen
        fprintf('\n  %2d   %3d   %2d \n', [it, res])
        
        vec_res(it) = res;
        if it > 1
        if (vec_res(it) > vec_res(it-1) && vec_res(it) < 1e-7)
            break
        end
        end
    end
    fprintf('\n Total iterations: %d \n\n', it)
    
    final_err = res; 
end