function omega = sabino_approx(m, M, ell)

s = (m/M);
kk = 1-(m/M)^2;
% K = 1.4 - 2.3*log10(s); 
K = log(4) - 0.5*log(s^2);
t = (2*(1:ell)-1)*K/(2*ell); %keyboard
[sn,cn,dn]=ellipj(t,kk); 
omega =sort(dn*M)';

end