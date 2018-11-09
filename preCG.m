function [xdisp, k] = preCG(x0, A, b, C, TOL)
%x0 & b need to be column vectors
%A needs to be symmetric
%C is the preconditioner s.t. M = C^2
k = 0;
x = x0;
r = b;
r0 = b;
z0 = C\r;
while norm(r, inf)>TOL
    k = k+1;
    z = C\r;
    if k == 1 
        p = z;
    else
        beta = (z'*r)/(z0'*r0);
        p = z + beta*p;
    end
    r0 = r;
    z0 = z;
    alpha = (z'*r)/(p'*A*p);
    x = x + alpha*p;
    r = r - alpha*A*p;       
end
xdisp = x;
end