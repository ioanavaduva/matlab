n = 2.^(2:10);
poles = zeros(9, 1);

for i = 1:9
    
    h = 1/n(i); 
    A = (diag(2*ones(n(i), 1)) + diag(-1*ones(n(i)-1, 1), 1) + diag(-1*ones(n(i)-1, 1), -1))/h^2;
    
    rhs1 = ones(n(i), 1);
    
    opts.tol=1e-4;
    emin = eigs(A, 1, 'smallestabs', opts);
    emax = eigs(A, 1,'largestabs',opts);
    
    m = 16; % number of poles; can change
    xi = emin + (emax-emin)*rand(1,m);
    [shifts, its, pole_norm] = irka_shifts(A, rhs1, xi, 1e-2);
    
    max_pole = max(shifts);
    
    poles(i) = max_pole;   
    
end

plot(log2(n), log(poles), 'x');