function [xdisp, k] = transfCG(x0, A, b, TOL)
%x0 & b need to be column vectors
%A needs to be symmetric
k = 0;
x = x0;
r0 = b;
r = b;
while norm(r, inf)>TOL
    %r ~= zeros(length(b), 1) not working - may be a bit too precise
    %norm(r, inf)>0 not working - I think it eventually divides by 0 and
    %gives all NaNs
    %>TOL makes sense as we want to get as close to 0 as we can but not
    %really reach it
    k = k+1;
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
end
xdisp = x;
disp(k);
end


