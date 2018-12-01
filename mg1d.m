% Multigrid method for 1D Poisson problem with damped Jacobi pre and post
% smoothing

function x = mg1d(n, b, T, w, maxit)
    % n is the number of unknowns
    % T is matrix of coefficients and b is rhs vector
    % w is the damping coefficient for damped_jacobi
    % maxit is the number of iterations to be used in damped_jacobi
    
    % pre-smoothing with damped Jacobi (do 50 iterations only)
    x = damped_jacobiM(w, zeros(length(b), 1), T, b, 10^-7, maxit);

    % compute the residual
    res = b - T*x;

    % generate restriction matrix
   
    k = log2(n+1);

    N = 2^(k-1)-1;

    RE = zeros(N,n);

    for i = 1:N
       RE(i,2*i-1:2*i+1) = [1 2 1]; 
    end

    RE = RE/4;

    % generate interpolation matrix
    II = 2*RE';

    % transfer residual to coarse grid; v is coarse grid residual
    v = zeros(N, 1);

    for i = 1:N
        v(i) = (res(2*i-1) + 2*res(2*i) + res(2*i+1))/4;
    end

    % transfer matrix T to coarse grid
    TC = RE * T * II;

    % solve residual equation to find error
    err = TC\v;

    % transfer error to fine grid; r is fine grid error
    erf = zeros(length(b), 1);

    erf(1) = err(1)/2;

    for j = 1:n/2
        erf(2*j) = err(j);
    end

    for j = 1:n/2-1
        erf(2*j+1) = (err(j) + err(j+1))/2;
    end

    erf(n) = err(length(err));

    % correct approximation (initial guess for damped Jacobi)
    x = x + erf;

    % post-smoothing Jacobi (50 iterations)
    x = damped_jacobiM(w, x, T, b, 10^-7, maxit);
end

