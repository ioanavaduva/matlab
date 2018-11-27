function x = multigrid(w, x0, A, b, TOL, maxit)
% maxit should be small as it is only needed for damped Jacobi and we want
% small number of iterations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate Restriction matrix R
    rep_vec = [1 2 1];
            
    R = zeros()
    R = (1/4)*
% generate Interpolation matrix I
    I = 2*R';
% pre-smoothing with damped Jacobi    
    x = damped_jacobiM(w, x0, A, b, TOL, maxit);
 % compute residual    
    res = b - A*x;
 % transfer residual to coarse grid
    res = R*res;
 % solve residual equation for error
    err = A\res;
 % transfer error to fine grid
    err = I*err;
 % correct approximation
    x = x + err;
 % post-smoothing with damped Jacobi
    x = damped_Jacobi(w, x, A, b, TOL, maxit);            
end
