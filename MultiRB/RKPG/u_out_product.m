function u = u_out_product(z, svec,iter,emin,emax)
% u_out computes u_A from equation 2.10 in [1]
% [1] BECKERMANN, BERNHARD. �AN ERROR ANALYSIS FOR RATIONAL GALERKIN 
% PROJECTION APPLIED TO THE SYLVESTER EQUATION.� SIAM Journal on Numerical
% Analysis, vol. 49, no. 5/6, 2011, pp. 2430�2450. 

% 
%     setup_Poisson_rank1rhs; % set-up of Poisson/convection-diffusion problem

    %opts.tol=1e-4;
    %emin = 1e-6; % min eigenvalue
    %emax=eigs(M{2},M{1},1,'LA',opts); % max eigenvalue

%     s = logspace(log10(emin)+0.1,log10(emax)-0.1, 6); % might need to change to take in the Zolotarev poles 

    u = ones(size(z));
    npoles = length(svec);
    svec = -svec;
    
%     display(svec)
    for k = 1:iter
        index = mod(k,npoles);
        if index == 0, index = npoles; end
        %fprintf('k = %i, index = %i\n',k,index);
        s = svec(index);
%     u_inner = abs((sqrt((z - emax)./(z - emin))*sqrt((s(6) - emax)./(s(6) - emin)) - 1)./((sqrt((z - emax)./(z - emin))*sqrt((conj(s(6)) - emax)./(conj(s(6)) - emin)) + 1)));
%     u = u.*u_inner;
% disp((s - emin)./(s - emax))
% disp(sqrt((s - emin)./(s - emax)))
% numerator = (sqrt((z - emax)./(z - emin))*sqrt((s - emin)./(s - emax)) - 1);
% denominator = (sqrt((z - emax)./(z - emin))*sqrt((conj(s) - emin)./(conj(s) - emax)) + 1);
% fprintf('numerator = %g, denominator = %g\n',numerator,denominator);
    
        u_inner = abs((sqrt((z - emax)./(z - emin))*sqrt((s - emin)./(s - emax)) - 1)./(sqrt((z - emax)./(z - emin))*sqrt((conj(s) - emin)./(conj(s) - emax)) + 1));
        u = u.*u_inner;
    end
    
    u = -u;
    
end
 
