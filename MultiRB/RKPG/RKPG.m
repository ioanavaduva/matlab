 function [X1, X2, final_err, vec_res, it, inner_it, avg_inner, error_vec] = RKPG(A, rhs1, rhs2, poles, tol, maxit, Xex_mat)
% Rational Krylov Subspace solver using the Petrov-Galerkin orthogonality
% condition. We currently solve Lyapunov equation XA + AX = rhs1*rhs2 with
% plan to extend to convection-diffusion matrix equation.
% !!! for Beckermann bound need to add extra output 'upper_vec'. 
% !!! for bounds on the error need extra input 'Xex_mat' and output
% 'error_vec'
% !!! for Ritz values need e_Ap in output which collects eigenvalues of
% projected matrix Ap

    inner_it = 0;
    i = 0; 
    res = 1; 
    it = 0; % the iterations also act as space dimension
    v0 = rhs1/norm(rhs1);
    V = v0; 
    vec_res(1) = res;
    %%%%  for Error plots error initialization
    error = 1;
    error_vec(1) = error;
    %%%%

%     %%%% for Beckerman bound --- can comment out when not interested in it
%     opts.tol=1e-4;
%     emin2= 1e-6; % min eigenvalue
%     emax2=eigs(A, 1,'LA',opts);
%     %%%%
    
    fprintf(' no.its  residual   no.inner its \n')
    
    while (res > tol && it < maxit)
        it = it + 1;
        i = i+1; 
        if i > length(poles) % cycle through poles          
            i = 1; 
        end
        
%         %%%%% compute the Beckermann upper bound -- very costly --- can comment
        % out when not interested in plotting
%         [val,fval,exitflag] =  fminbnd(@(z) u_out_product(z, poles, it, emin2, emax2), -emax2, -emin2);
%         fval = -fval;
%         const = 4 + 4*sqrt(2*cond(A));
%         upper_bound = const*fval;
%         upper_vec(it+1) = upper_bound;
%         %%%%%
        
        % choose basis 
        V = get_rk_basis(A, poles(i), V); %keyboard % generate the rational Krylov basis
%         V = get_poly_basis(A, V); % generate the polynomial (standard) Krylov basis
%         V = get_ek_basis(A, V);
        
        % project matrix A and rhs1/2
        Ap = V'*A*V; 
        rhs1p = V'*rhs1;
        rhs2p = V'*rhs2;
%         e_Ap{it} = eig(Ap); (uncomment for Ritz values)
        
        % solve projected problem
         Y = lyap(-Ap, rhs1p*rhs2p'); 
%         proj_dim = size(A, 2);
%         tol_inner = tol*1e-1;
%         y0 = zeros(proj_dim, proj_dim);
%         App = mat2cell(Ap, [size(Ap, 2), 0]);
%         [Y, inner_it] = cgkron(App, App, rhs1p*rhs2p', y0, proj_dim*proj_dim, tol_inner);

%         norm(Ap*Y+Y*Ap'-rhs1p*rhs2p', 'fro')/norm(rhs1p*rhs2p', 'fro')

        % obtain low-rank factors X_1 & X_2 using SVD
        [uu, ss, vv] = svd(Y, 0);
        X1 = V*uu*ss;
        X2 = vv'*V';
        
        %%% Exact solution at each iteration (only need for bounds & to compare 1-&2-sided proj)
        XX = X1*X2; %keyboard
        error = norm(Xex_mat - XX);
        error_vec(it) = error;
        %%%
        
        % project back
%         X_hat = V*Y*V';
        
        % compute residual
        res = norm(X1*X2*A + A*X1*X2 - rhs1*rhs2', 'fro')/norm(rhs1*rhs2', 'fro');
        
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

