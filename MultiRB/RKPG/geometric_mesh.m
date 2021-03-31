function [x,dx] = geometric_mesh(n,r)

if r == 1
    dx = 1/(n+1)*ones(1,n);
    x = 0:dx:1;
else
    
    a = 0.5*(1-r)/(1-r^((n-2)/2+1));
    
    temp = a*r.^(0:(n-2)/2);
%     keyboard;
    dx = [temp fliplr(temp)];
    x = [0 cumsum(dx)];
    
end
% end

