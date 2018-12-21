% Multigrid method for 1D Poisson problem with damped Jacobi pre and post
% smoothing

function x = vcycle1d(n, x, b, T, w, levels, resC ,Ts, xst)
    
L = levels;
Ts{L} = T;
resC{L} = b;
if levels ~= 1
% resC = b; Ts = T, xst = zeros (just like x)
    % n is the number of unknowns of the form n = 2^k-1
    % T is matrix of coefficients and b is rhs vector
    % w is the damping coefficient for damped_jacobi
    % maxit is the number of iterations to be used in damped_jacobi
    % pre-smoothing with damped Jacobi (do 3 iterations only)
    
    
    x = damped_jacobiM(w, x, T, b, 10^-7, 3);
    xst{L} = x;
    % compute the residual
    res = b - T*x;
    
 %   if norm(res) > TOL
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
        resC{L-1} = v;
        
                
        % transfer matrix T to coarse grid
        TC = RE * T * II;
        Ts{L-1} = TC;
        n = N;
        
        err = vcycle1d(N, zeros(length(v), 1), v, TC, w, levels-1, resC, Ts, xst);
        else
        % solve residual equation to find error
            err = Ts{L}\resC{L};
        end
        
        
        for L = 1:levels-1

                N = n;
                
                k = log2(n+1);
                n = 2^(k+1)-1;
                
                
                % transfer error to fine grid; erf is fine grid error
                erf = zeros(length(resC{L}), 1);        
                erf(1) = err(1)/2;        
                for j = 1:n/2
                    erf(2*j) = err(j);
                end        
                for j = 1:n/2-1
                    erf(2*j+1) = (err(j) + err(j+1))/2;
                end        
                erf(n) = err(length(err));
                size(erf)
                size(xst{L+1})
                % correct approximation (initial guess for damped Jacobi in post-smoothing)
                x = erf + xst{L+1};
                
                size(Ts{L+1})
                size(resC{L+1})
                size(x)
                %post-smoothing Jacobi (3 iterations)
                x = damped_jacobiM(w, x, Ts{L+1}, resC{L+1}, 10^-7, 3);

                err = x;                 
            end
%        % transfer error to fine grid; r is fine grid error
%        erf = zeros(length(b), 1);
        
%       erf(1) = err(1)/2;

%        for j = 1:n/2
%            erf(2*j) = err(j);
%        end
        
%        for j = 1:n/2-1
%            erf(2*j+1) = (err(j) + err(j+1))/2;
%        end
        
%        erf(n) = err(length(err));
        
        % correct approximation (initial guess for damped Jacobi)
%        x = x + erf;
        
        % post-smoothing Jacobi (3 iterations)
%        x = damped_jacobiM(w, x, T, b, 10^-7, 3);

end



