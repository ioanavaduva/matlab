
addpath(genpath('../../rktoolbox'));

vec_n = [50, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400];
resid_mat = zeros(101, length(vec_n));
tol = 1e-14;
maxit = 100;

for i = 1:length(vec_n)
    n = vec_n(i);
    
    % create matrix & rhs
    h = 1/n; eps = 1;
    A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;
    rhs1 = ones(n, 1);
    rhs2 = ones(n, 1);
    
    % Get smallest and largest eigenvalues
    emin = 1e-6; 
    opts.tol=1e-4;
    emax = eigs(A, 1,'LA',opts);
    
    % 4 positive imaginary parts of Zolotarev poles
    bb = emax - emin + 1;

    k = 4;      % rational degree
    b = bb;     % sign function on [-10,-1]\cup [1,10]
    r = rkfun.gallery('sign', k, b);
    % poles(r)
    po = imag(poles(r));
    poles_Zolo = po(po >= 0);
    
    % Solve with RKPG.m
    [~, ~, ~, vec_res, ~, ~, ~] = RKPG(A, rhs1, rhs2, poles_Zolo, tol,  maxit);
    
    resid_mat(:,i) = vec_res';
    
    i = i+1;
end
    
