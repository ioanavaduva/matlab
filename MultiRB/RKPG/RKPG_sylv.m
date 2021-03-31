 function [X1, X2, final_err, vec_res, it, inner_it, avg_inner] = RKPG_sylv(A, T, rhs1, rhs2, poles, tol, maxit)

    inner_it = 0;
    i = 0; 
    res = 1; 
    it = 0; % the iterations also act as space dimension
    v0 = rhs1/norm(rhs1);
    V = v0; 
    vec_res(1) = res;
    
    fprintf(' no.its  residual   no.inner its \n')
    
    while (res > tol && it < maxit)
        it = it + 1;
        i = i+1; 
        if i > length(poles) % cycle through poles          
            i = 1; 
        end
        
        % choose basis 
        V = get_rk_basis(A, poles(i), V); 
        
        % project matrix A and rhs1/2
        Ap = V'*A*V; 
        rhs1p = V'*rhs1;
        
        % solve projected problem
         Y = lyap(-Ap, -T,  rhs1p*rhs2'); 

        % obtain low-rank factors X_1 & X_2 using SVD
        [uu, ss, vv] = svd(Y, 0);
        X1 = V*uu*ss;
        X2 = vv';
        
        % compute residual
        res = norm(X1*X2*T + A*X1*X2 - rhs1*rhs2', 'fro')/norm(rhs1*rhs2', 'fro');
        
        % Print details to screen
        fprintf('\n  %2d   %3d   %2d \n', [it, res, inner_it])
        
        vec_res(it) = res;
        if it > 1
        if (vec_res(it) > vec_res(it-1) && vec_res(it) < 1e-8)
            break
        end
        end
    end
    fprintf('\n Total iterations: %d \n\n', it)
    
    final_err = res;
    avg_inner = mean(inner_it);
end


