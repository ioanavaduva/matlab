% Multigrid method for 1D Poisson problem with damped Jacobi pre and post
% smoothing; V-cycle with j levels (usually 4)

function [x, it] = Vmg1d(n, b, T, w, maxit, TOL)
    % n is the number of unknowns of the form n = 2^k-1
    % T is matrix of coefficients and b is rhs vector
    % w is the damping coefficient for damped_jacobi
    % maxit is the number of iterations to be used in damped_jacobi
    
    it = 0; 
    TC = T;
    resC = b;
    N = n;
    
    for it = 1:maxit
        
        for L = 1:1:2
                        
            x0 = zeros(length(resC),1);
            
            % pre-smoothing with damped Jacobi (do 3 iterations only)
            x = damped_jacobiM(w, x0, TC, resC, 10^-7, 3);

            % compute the residual
            res = resC - TC*x;
            
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
            
            % transfer matrix T to coarse grid
            TC = RE * TC * II;
            
            % transfer residual to coarse grid; resC is coarse grid residual
            
            resC = zeros(N, 1);
            
            for i = 1:N
                resC(i) = (res(2*i-1) + 2*res(2*i) + res(2*i+1))/4;
            end
            
            n = N;
            L = L+1;          
        end
            resC  
            n
            TC
            N
            x
        
        % solve residual equation to find error
        err = TC\resC
            
        for L = 2:-1:1
            
            x0 = zeros(length(resC),1); 
            
            N = n
            
            k = log2(n+1)
            
            n = 2^(k+1)-1
                       
            % transfer error to fine grid; erf is fine grid error
            erf = zeros(length(resC), 1);
            
            erf(1) = err(1)/2;
            
            for j = 1:n/2
                erf(2*j) = err(j);
            end
            
            for j = 1:n/2-1
                erf(2*j+1) = (err(j) + err(j+1))/2;
            end
            
            erf(n) = err(length(err));
            
            erf           

            % correct approximation (initial guess for damped Jacobi)
            
            x = x + erf         
            
            % generate restriction matrix            
            
            RE = zeros(N,n);
            
            for i = 1:N
                RE(i,2*i-1:2*i+1) = [1 2 1]; 
            end
            
            RE = RE/4;
            
            % generate interpolation matrix
            
            II = 2*RE';
            
            % transfer matrix T to fine grid
            TC = II * TC * RE;
            
            % post-smoothing Jacobi (3 iterations)
            x = damped_jacobiM(w, x, TC, resC, 10^-7, 3)
           
            err = x;
            L = L + 1;                
        end
        it = it+1; 
    end
end