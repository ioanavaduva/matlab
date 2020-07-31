function [X1, X2, vec_res, it, final_err, upper_vec] = RKPGblock2(A, rhs1, rhs2,poles, tol, maxit)
% For Beckermann bound need extra ''upper_vec''
    i = 1; j=i;
    res = 1; 
    it = 1; 
    vec_res(1) = res;
    V = get_start_basis(rhs1); % computes the starting basis block of same dimension as rhs)
    
    %%%% for Beckerman bound --- can comment out when not interested in it
    opts.tol=1e-4;
    emin2= 1e-6; % min eigenvalue
    emax2=eigs(A, 1,'LA',opts);
    %%%%
    
    fprintf(' no.its  residual   basis size \n')
    while(res > tol && it < maxit)
        if i > length(poles)
            i=1;
        end
        for k = 1:length(rhs1(1, :)) % k refers to what column we're filling in
            col = (j-1)*length(rhs1(1, :)) + k; %keyboard
            V = get_rk_block_basis(A, poles(i), V, col);
        end
            % project matrix A and rhs1/2
            Ap = V'*A*V; 
            rhs1p = V'*rhs1;
            rhs2p = V'*rhs2;
            % solve projected problem
            Y = lyap(-Ap, rhs1p*rhs2p'); 
            % obtain low-rank factors X_1 & X_2 using SVD
            [uu, ss, vv] = svd(Y, 0);
            X1 = V*uu*ss;
            X2 = vv'*V';
            % compute approx solution
            XX = X1*X2;
            % compute residual
            %%%%% compute the Beckermann upper bound -- very costly --- can comment
            % out when not interested in plotting
            [val,fval,exitflag] =  fminbnd(@(z) u_out_product(z, poles, it, emin2, emax2), -emax2, -emin2);
            fval = -fval;
            const = 4 + 4*sqrt(2*cond(A));
            upper_bound = const*fval;
            upper_vec(it) = upper_bound;
            %%%%%
            res = norm(XX*A + A*XX - rhs1*rhs2', 'fro')/norm(rhs1*rhs2', 'fro');
            % Print details to screen
            fprintf('\n  %2d   %3d   %2d \n', [it, res, length(V(1, :))])
            vec_res(it) = res;
            if res < tol || it > maxit
                break;
            end
            it = it+1;
        i = i+1; j=j+1;
    end
    final_err = res;
end