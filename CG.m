function [x, iter] = CG(x0, A, b, TOL, maxit)
%x0 & b need to be column vectors
%A needs to be symmetric
x = x0;
r0 = b;
r = b;
iter = 0;
for k = 1:maxit
    if norm(r, inf) > TOL
        if k == 1
            p = r;
        else
            beta = (r'*r)/(r0'*r0);
            p = r + beta*p;
        end
        r0 = r;
        alpha = (r'*r)/(p'*A*p);
        x = x + alpha*p;
        r = r - alpha*A*p;       
    else
        break;
    end
    iter = iter + 1;
end
disp(iter);
end


