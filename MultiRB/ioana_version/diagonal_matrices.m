% Diagonal matrices to look at eigenvalues 

%%%%%%%%%%% SETUP

n = 100; 

% Uniform spacing
x1 = linspace(1e-8, n^2, n);
D1 = diag(x1);

% Logspace 
x2 = logspace(1e-6, 5, n);
D2 = diag(x2);

% Clustered in centre
x3 = [linspace(1e-8, n/4, n/4), linspace(n/3, 2*n/3, 3*n/4-n/4), linspace(3*n/4, n, n-3*n/4)];
D3 = diag(x3);

% Smallest eigenvalue away from the rest
x4 = [1e-8, logspace(3, 5, n-1)];
D4 = diag(x4);

A = D4; I = speye(n);
M0 = I; N0 = A; 

M = {M0, N0};
N = {N0, M0};
% rank-one symmetric rhs
rhs1 = ones(n,1);
rhs2=rhs1; 
% rhs=ones(n^2,1);

% rank-one nonsymmetric rhs
% rhs1 = ones(n, 1);
% rhs2 = -rhs1 + [0, 0, 1, 0, 0]';
% NNn = 0.0001;
% rhs2 = NNn*randn(n, 1) + ones(n, 1)
% rhs = rhs1*rhs2'; rhss = rhs(:);

% rhs1 = cos(2*n*pi*x)';
% rhs2 = sin(2*n*pi*x)';
rhs = rhs1*rhs2';
rhss = rhs(:);

% preconditioner
P1 = speye(n); 
P = P1*P1';

% kronecker product form of matrix eq
AA  = kron(A', I) + kron(I', A); % + kron(N1', M1) + kron(N2', M2);

%%%%%%%%%%%%%% DRIVER

fprintf('\n Dim_M = %d,  Dim_N = %d, no. terms = %d  \n', size(M{1},1),size(N{1},1),length(M))
fprintf('\n No. equations (Ax=b): %d  \n', size(M{1},1)*size(N{1},1))

n=size(M{1},1); m=size(N{1},1);

% ----------------Choose shifts ------------------------------------
opts.tol=1e-4;  % NB: 'eigs' is VERY sensitive to this.

% This could be done more cheaply with data from a COARSER problem
% min and max eigs of hat{K}_{1} - 1st term could be strictly positive
emin2=eigs(M{2},M{1},1,'SA',opts);
emax2=eigs(M{2},M{1},1,'LA',opts);
n_m=size(M,2);
emean=mean([(emin2),abs(emax2)]);
alphas(1)=1-emean;
alphas(2:n_m)=(1)*ones(n_m-1,1);

% Shifts given by 2 + 2cos(j*pi/(k+1))
% for j = 1:6
%     shifts(j) = eps*(2 + 2*cos((j-1)*pi/6))/h^2;
% end
% s_parameter = shifts;

% Shifts that are the exact eigenvalues
% mm = round(linspace(1, n, 6));
% ei = sort(eig(A));
% for i = 1:6
%     shifts(i) = ei(mm(i));   
% end
% s_parameter = shifts;

% 'parameter-free version' for above alphas;
%     s_parameter(1)=2-alphas(1);
%     s_parameter(2)=1;
%     s_parameter                             % Uncomment to print to screen

% multiple parameter strategy - [lambda_min (M), lambda_max(M), k); using
% elliptic functions

aa=emin2;%+alphas(1); 
bb=emax2;%+alphas(1);
S_interval=[aa,bb];
s_nodes = 6;                           % Choose 2 nodes (could vary)
snew = get_nodes2(aa,bb,s_nodes);      % Use interval for A_1;
s_parameter=snew;

% ------------ Options for MultiRB --------------------------------
res_method='4';                  % can't change this
rat_solve='1'; mmax=200;         % can change
param.max_space_dim=mmax;
param.period=1;                  % can change
param.rat_solve=rat_solve;
param.res_method=res_method;
tic;
fprintf('\n -------------------- MultiRB Solve ------------------------\n \n')

[X1,X2,dimV,final_err,avg_inner,error_vec,iv_vec]=MultiRB_Poisson_rank1rhs_2sided(M,N,rhs1,rhs2,P,P1,param,s_parameter);
%X1(psort,:)=X1;
etoc=toc; 
fprintf('\n Total execution time: %9.4e seconds \n',etoc)
fprintf('\n ----------------------------------------------------------\n \n')
fprintf('no_terms   dim_M   dim_N   n_k   Rank    final_err   avg_inner  time(s)  \n')
fprintf('\n  %2d       %3d   %d    %3d    %2d    %9.4e    %4.2f     %9.4e  \n \n', [2, n, m, dimV, size(X1,2), final_err, avg_inner, etoc])
fprintf('\n----------------------------------------------------------\n \n')
