function [x, iter] = steepd(x0, A, b, TOL, maxit)
%x0 needs to be column vector
%b needs to be row vector
%A needs to be symmetric
iter = 0;
x = x0;
r = zeros(length(b),1);
%xold = x0;
%r = b';
for k = 1:maxit
     r = b - A * x;
     if norm(r, inf) > TOL %want to minimise the norm of r b/c if r is close 
         %to 0 then b'-AX is close to being 0 which satisfies Ax=b (our
         %system that we want to solve)
        a = (r'*r)/(r'*A*r);
        x = x0 + a * r;  
     else 
         break;
     end
    x0 = x;
    iter = iter+1;
end
end
