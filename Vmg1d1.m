% Multigrid method for 1D Poisson problem with damped Jacobi pre and post
% smoothing

function x1 = Vmg1d1(n, b, T, w)
    % n is the number of unknowns of the form n = 2^k-1
    % T is matrix of coefficients and b is rhs vector
    % w is the damping coefficient for damped_jacobi
    % maxit is the number of iterations to be used in damped_jacobi
  
    x = zeros(length(b),1);
    
    % pre-smoothing with damped Jacobi (do 3 iterations only)
    x1 = damped_jacobiM(w, x, T, b, 10^-7, 3);
    
    % compute the residual
    res = b - T*x1;
 
    % generate restriction matrix
    k = log2(n+1);
    N = 2^(k-1)-1;
    
    RE1 = zeros(N,n);
    for i = 1:N
        RE1(i,2*i-1:2*i+1) = [1 2 1]; 
    end
    RE1 = RE1/4;
    
    % generate interpolation matrix
    II1 = 2*RE1';
    
    % transfer residual to coarse grid; v is coarse grid residual
    v1 = zeros(N, 1);
    for i = 1:N
        v1(i) = (res(2*i-1) + 2*res(2*i) + res(2*i+1))/4;
    end
    
    % transfer matrix T to coarse grid
    TC = RE1 * T * II1;
    
    % repeat
    x = zeros(length(v1),1);
    
    % pre-smoothing with damped Jacobi (do 3 iterations only)
    x2 = damped_jacobiM(w, x, TC, v1, 10^-7, 3);
    
    % compute the residual
    res = v1 - TC*x2;
    
    n = N; 
    % generate restriction matrix
    k = log2(n+1);
    N = 2^(k-1)-1;
    RE2 = zeros(N,n);
    for i = 1:N
        RE2(i,2*i-1:2*i+1) = [1 2 1]; 
    end
    RE2 = RE2/4;
    
    % generate interpolation matrix
    II2 = 2*RE2';
    
    % transfer residual to coarse grid; v is coarse grid residual
    v2 = zeros(N, 1);
    for i = 1:N
        v2(i) = (res(2*i-1) + 2*res(2*i) + res(2*i+1))/4;
    end
    
    % transfer matrix T to coarse grid
    TC2 = RE2 * TC * II2;
    
    % repeat
    x = zeros(length(v2),1);
    
    % pre-smoothing with damped Jacobi (do 3 iterations only)
    x3 = damped_jacobiM(w, x, TC2, v2, 10^-7, 3);
    
    % compute the residual
    res = v2 - TC2*x3;
    
    n = N; 
    % generate restriction matrix
    k = log2(n+1);
    N = 2^(k-1)-1;
    RE3 = zeros(N,n);
    for i = 1:N
        RE3(i,2*i-1:2*i+1) = [1 2 1]; 
    end
    RE3 = RE3/4;
    
    % generate interpolation matrix
    II3 = 2*RE3';
    
    % transfer residual to coarse grid; v is coarse grid residual
    v3 = zeros(N, 1);
    for i = 1:N
        v3(i) = (res(2*i-1) + 2*res(2*i) + res(2*i+1))/4;
    end
    
    % transfer matrix T to coarse grid
    TC3 = RE3 * TC2 * II3;
    
    % solve residual equation to find error
    err = TC3\v3;
    
    % transfer error to fine grid; erf is fine grid error
    erf = zeros(length(v2), 1);
    erf(1) = err(1)/2;
    for j = 1:n/2
        erf(2*j) = err(j);
    end
    for j = 1:n/2-1
        erf(2*j+1) = (err(j) + err(j+1))/2;
    end
    erf(n) = err(length(err));
    
    % correct approximation (initial guess for damped Jacobi)
    x3 = x3 + erf;
    
    % post-smoothing Jacobi (3 iterations)
    x3 = damped_jacobiM(w, x3, TC2, v2, 10^-7, 3);
    
    % transfer error to fine grid; erf is fine grid error
    erf = zeros(length(v1), 1);
    erf(1) = err(1)/2;
    for j = 1:n/2
        erf(2*j) = err(j);
    end
    for j = 1:n/2-1
        erf(2*j+1) = (err(j) + err(j+1))/2;
    end
    erf(n) = err(length(err));
    
    % correct approximation (initial guess for damped Jacobi)
    x2 = x2 + erf;
    
     % post-smoothing Jacobi (3 iterations)
    x3 = damped_jacobiM(w, x2, TC, v1, 10^-7, 3);
    
    % transfer error to fine grid; erf is fine grid error
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
    x1 = x1 + erf;
    
     % post-smoothing Jacobi (3 iterations)
    x1 = damped_jacobiM(w, x1, T, b, 10^-7, 3);
    
end
