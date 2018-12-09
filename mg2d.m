% Multigrid method for 2dD Poisson problem with damped Jacobi pre and post
% smoothing

function [x, it] = mg2d(n, b, A, w, maxit, TOL)
    % n is the number of unknowns of the form n = 2^k-1
    % T is matrix of coefficients and b is rhs vector
    % w is the damping coefficient for damped_jacobi
    % maxit is the number of iterations to be used in damped_jacobi
    it = 0;
    x = zeros(length(b),1);
    for it = 1:maxit
        % pre-smoothing with damped Jacobi (do 3 iterations only)
        x = damped_jacobiM(w, x, A, b, 10^-7, 3);
        
        % compute the residual
        res = b - A*x;
        
        if norm(res) > TOL
            % generate restriction & interpolation matrices

            k = log2(n+1);

            N = 2^(k-1)-1;

            RE = zeros(N,n);

            for i = 1:N
               RE(i,2*i-1:2*i+1) = [1 2 1]; 
            end

            RE = RE/4;
            II = 2*RE';
            
            II2d = kron(II, II);
            
            RE2d = kron(RE, RE);
            
            % make vector of residuals into matrix with N columns
            res = reshape(res, [n, n]);
            % transfer residual to coarse grid; v is coarse grid residual
            for i = 1:(n+1)/2-1
                for j = 1:(n+1)/2-1
                    v(i, j) = (res(2*i-1, 2*j-1) + res(2*i-1, 2*j+1) + res(2*i+1, 2*j-1) + res(2*i+1, 2*j+1) + 2*(res(2*i, 2*j-1)+res(2*i, 2*j+1) + res(2*i-1, 2*j) + res(2*i+1, 2*j)) + 4*res(2*i, 2*j))\16;
                end
            end
            
            v = reshape(v, [N, 1]);

            % transfer matrix T to coarse grid
            AC = RE2d*A*II2d;

            % solve residual equation to find error
            err = AC\v;

            % transfer error to fine grid; r is fine grid error
            erf = zeros(length(b), 1);

            erf(1,:) = err(1)/2;
            for i = 1:n
                erf(2*i, :) = err(i, :);
                erf(2*i+1, :) = (err(i, :)+ err(i+1, :))/2;
            end
            
            erf(:, 1) = erf(:, 1)/2;
            for j = 1:n
                erf(:, 2*j) = erf(:, j);
                erf(:, 2*j+1) = (erf(:, j)+ erf(:, j+1))/2;
            end

            % correct approximation (initial guess for damped Jacobi)
            x = x + erf(:, 1);

            % post-smoothing Jacobi (3 iterations)
            x = damped_jacobiM(w, x, A, b, 10^-7, 3);
        
        else
            break;
        end
        it=it+1;     
    end
end
