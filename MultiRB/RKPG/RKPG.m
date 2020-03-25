function [X1, X2, final_err, vec_res, it, inner_it, avg_inner] = RKPG(A, rhs1, rhs2, poles, tol, maxit)
% Rational Krylov Subspace solver using the Petrov-Galerkin orthogonality
% condition. We currently solve Lyapunov equation XA + AX = rhs1*rhs2 with
% plan to extend to convection-diffusion matrix equation.

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
        if i > length(poles) % cycle through poles                                                              v
            i = 1; 
        end
        V = get_rk_basis(A, poles(i), V); %keyboard % generate the rational Krylov basis
%         rank(V)
%         V = get_poly_basis(A, V); % generate the polynomial (standard) Krylov basis
%         V = get_ek_basis(A, V);
        
        % project matrix A and rhs1/2
        Ap = V'*A*V; %keyboard
        rhs1p = V'*rhs1;
        rhs2p = V'*rhs2;
        
        % solve projected problem
        
         Y = lyap(-Ap, rhs1p*rhs2p'); %keyboard
%         proj_dim = size(A, 2);
%         tol_inner = tol*1e-1;
%         y0 = zeros(proj_dim, proj_dim);
%         App = mat2cell(Ap, [size(Ap, 2), 0]);
%         [Y, inner_it] = cgkron(App, App, rhs1p*rhs2p', y0, proj_dim*proj_dim, tol_inner);

%         norm(Ap*Y+Y*Ap'-rhs1p*rhs2p', 'fro')/norm(rhs1p*rhs2p', 'fro')

        % obtain low-rank factors Y_1 & Y_2 using SVD
        [uu, ss, vv] = svd(Y, 0);
        X1 = V*uu*ss;
        X2 = vv'*V';
        
        % project back
%         X_hat = V*Y*V';
        
        % compute residual
        res = norm(X1*X2*A + A*X1*X2 - rhs1*rhs2', 'fro')/norm(rhs1*rhs2', 'fro');
        
        % Print details to screen
        fprintf('\n  %2d   %3d   %2d \n', [it, res, inner_it])
        
        vec_res(it + 1) = res;
    end
    
    fprintf('\n Total iterations: %d \n\n', it)
    
    final_err = res;
    avg_inner = mean(inner_it);
end

