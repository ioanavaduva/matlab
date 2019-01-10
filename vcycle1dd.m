function [x, Out] = vcycle1dd(n, x, b, bb, T, TT, w, levels, lvl, Out)  
    
    x = damped_jacobiM(w, x, T, b, 1e-7, 3);
    res = b - T*x;
    
    k = log2(n+1);
    N = 2^(k-1)-1;
    RE = zeros(N,n);
    q = [1,2,1] / 4;
    for i = 1:N
        RE(i, 2*i-1:2*i+1) = q;
    end
    
    II = 2 * RE';
    TC = RE * T * II;
    
    for i = 1:N   
        v(i) = (res(2*i-1) + 2*res(2*i) + res(2*i+1))/4;
    end
    
    Out(levels).v  = v';
    Out(levels).TC = TC;
    Out(levels).x = x;
  %  size(x)
    if levels ~= 1
        [x, Out] = vcycle1dd(N, zeros(length(v), 1), v', bb, TC, TT, w, levels-1, lvl, Out);
 %      a= size(err)
    else
        return;
    end
    Out(5).TC = TT;
    Out(5).v  = bb;

    er = Out(1).TC\Out(1).v;
    
    k = log2(n+1);
    N = 2^(k-1)-1;
    n = N;
    k = log2(n+1);
    N = 2^(k-1)-1;
    n = N;
    
    for L = 1:lvl
        
        N = n;
        k = log2(n+1);
        n = 2^(k+1)-1;
        % transfer error to fine grid; erf is fine grid error
        erf = zeros(length(er), 1);
        erf(1) = er(1)/2;
        for j = 1:n/2
            erf(2*j) = er(j);
        end
        for j = 1:n/2-1
            erf(2*j+1) = (er(j) + er(j+1))/2;
        end
        erf(n) = er(length(er));
        
        % correct approximation (initial guess for damped Jacobi in post-smoothing)
        x = erf + Out(L).x;
        
        %size(x)
        %size(Out(L+1).TC)
        %size(Out(L+1).v)
             
        %post-smoothing Jacobi (3 iterations)
        x = damped_jacobiM(w, x, Out(L+1).TC, Out(L+1).v, 10^-7, 3);
        er = x;
    end
end