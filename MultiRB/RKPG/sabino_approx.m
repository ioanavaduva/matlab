function omega = sabino_approx(m, M, ell)

s = (m/M);
kk = 1-(m/M)^2;
K = 1.4 - 2.3*log10(s); %keyboard %0.809079148555692*kk + 1.506421783135102; keyboard
t = (2*(1:ell)-1)*K/(2*ell); %keyboard
[sn,cn,dn]=ellipj(t,kk); 
omega =sort(dn*M)';

end