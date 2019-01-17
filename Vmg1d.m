% Nonrecursive V cycle multigrid method for 1D Poisson problem with damped 
% Jacobi pre and post smoothing

% n is the number of unknowns of the form n = 2^k-1
% T is matrix of coefficients and b is rhs vector
% w is the damping coefficient for damped_jacobi
% maxit is the number of iterations
% maxlev is the number of levels we go through
% TOL is tolerance level for convergenc: 1e-4 is ok here

function [x, iter] = Vmg1d(n, b, T, w, maxit, TOL, maxlev)

    T_original = T;
    resC = b;
    rst{1} = resC;
    iter = 0;
    x = zeros(length(b), 1);
    n_orig = n;
    for iter = 1:maxit
        TC = T_original;
        n = n_orig;
        r = b - T*x;
        x0 = x;
        resC = b;
        if norm(r) > TOL
            
            for L = 1:1:maxlev
                if L > 1
                x0 = zeros(length(resC),1);
                end
                % pre-smoothing with damped Jacobi (do 3 iterations only)
                x = damped_jacobiM(w, x0, TC, resC, 1e-7, 3);
                xst{L} = x;
                Tst{L} = TC;
                
                % compute the residual
                res{L} = resC - TC*x;
                
                    
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
                    resC(i) = (res{L}(2*i-1) + 2*res{L}(2*i) + res{L}(2*i+1))/4;
                end        
                rst{L+1} = resC;
                n = N;
            
            end
            
           
            % solve residual equation to find error
            err = TC\resC;
            
           
            for L = maxlev:-1:1

                N = n;
                k = log2(n+1);
                n = 2^(k+1)-1;
                

                % transfer error to fine grid; erf is fine grid error
                erf = zeros(n, 1); 
                erf(1) = err(1)/2;        
                for j = 1:n/2
                    erf(2*j) = err(j);
                end        
                for j = 1:n/2-1
                    erf(2*j+1) = (err(j) + err(j+1))/2;
                end        
                erf(n) = err(length(err))/2;
                


                % correct approximation (initial guess for damped Jacobi in post-smoothing)
                x = erf + xst{L};


                %post-smoothing Jacobi (3 iterations)
                x = damped_jacobiM(w, x, Tst{L}, rst{L}, 1e-7, 3);
                

                err = x;                 
            end
        else
            break;
        end
    end
end