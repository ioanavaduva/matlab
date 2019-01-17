% Multigrid method for 2D Poisson problem with damped Jacobi pre and post
% smoothing; V-cycle with j levels (usually 4)

function [x, iter] = Vmg2d(n, b, A, w, maxit, TOL, maxlev)
    % n is the number of unknowns of the form n = 2^k-1
    % T is matrix of coefficients and b is rhs vector
    % w is the damping coefficient for damped_jacobi
    % maxit is the number of iterations to be used in damped_jacobi

    A_original = A;
    resC = b;
    rst{1} = resC;
    iter = 0;
    x = zeros(length(b), 1);
    n_orig = n;
    for iter = 1:maxit
        AC = A_original;
        n = n_orig;
        r = b - A*x;
        x0 = x;
        resC = b;
        if norm(r) > TOL
            for L = 1:1:maxlev
                if L > 1
                x0 = zeros(length(resC),1);
                end
                % pre-smoothing with damped Jacobi (do 3 iterations only)
                x = damped_jacobiM(w, x0, AC, resC, 1e-7, 3);
             
                xst{L} = x;
                Ast{L} = AC;

                % compute the residual
                res{L} = resC - AC*x;
                
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
                
                
                II2d = kron(sparse(II), sparse(II));
            
                RE2d = kron(sparse(RE), sparse(RE));

                % transfer matrix T to coarse grid
                AC = sparse(RE2d) * sparse(AC) * sparse(II2d); 

                % transfer residual to coarse grid; resC is coarse grid residual
                res{L} = reshape(res{L}, [n,n]);
                resC = zeros(N,N);     
                for i = 1:(n+1)/2 - 1
                    for j = 1:(n+1)/2 - 1
                        resC(i, j) = (res{L}(2*i-1, 2*j-1) + res{L}(2*i-1, 2*j+1) + res{L}(2*i+1, 2*j-1) + res{L}(2*i+1, 2*j+1) + 2*(res{L}(2*i, 2*j-1)+res{L}(2*i, 2*j+1) + res{L}(2*i-1, 2*j) + res{L}(2*i+1, 2*j)) + 4*res{L}(2*i, 2*j))/16;
                    end
                end      
                resC = reshape(resC, [N^2,1]);
                
                rst{L+1} = resC;
                n = N;
            end
           
            % solve residual equation to find error
            err = AC\resC;
            
            for L = maxlev:-1:1

                N = n;
                k = log2(n+1);
                n = 2^(k+1)-1;

                % transfer error to fine grid; erf is fine grid error
                err = reshape(err, [N, N]);
                erf = zeros(n, n);
                
                for i = 1:n/2
                    for j = 1:n/2
                        erf(2*i, 2*j) = err(i, j);
                    end
                end
                for i = 1:n/2-1
                    for j = 1: n/2
                        erf(2*i+1, 2*j) = (err(i, j) + err(i+1, j))/2;
                    end
                end
                for i = 1:n/2
                    for j =1: n/2-1
                        erf(2*i, 2*j+1) = (err(i, j) + err(i, j+1))/2;
                    end
                end
                for i = 1:n/2-1
                    for j =1: n/2-1
                         erf(2*i+1, 2*j+1) = (err(i, j) + err(i+1, j) + err(i, j+1) + err(i+1, j+1))/4;
                    end
                end
                 
                erf = reshape(erf, [n^2, 1])
                
                % correct approximation (initial guess for damped Jacobi in post-smoothing)
                x = erf(:,1) + xst{L};

                %post-smoothing Jacobi (3 iterations)
                x = damped_jacobiM(w, x, Ast{L}, rst{L}, 1e-7, 3);

                err = x;                 
            end
        else
            break;
        end
    end
end