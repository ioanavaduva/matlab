function [x,h] = geometric_mesh3(n,r)

if mod(n,2)~=0
    error('n must be even')
end

h = 1/(n+1);
x = 0:h:1;

if r == 1
    h = diff(x);
else
    
    np = n/2+1;
    xl = x(n/2+1);
    xr = x(n/2+2);
    
    a = (1-xl)*(1-r)/(1-r^np);
    
    dx = cumsum(a*r.^(0:np-1));
    
    xpr = [xl + dx];
    xpl = [xr - dx(end:-1:1)];
    
    x = sort([xpl xpr]); x(1) = 0; x(end) = 1;
    h = diff(x);
end
