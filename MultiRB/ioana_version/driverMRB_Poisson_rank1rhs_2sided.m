% Driver for MultiRB.m code, version 1

clear all
addpath(genpath('../../rktoolbox'));

setup_Poisson_rank1rhs;
%setup_nonsym;

fprintf('\n Dim_M = %d,  Dim_N = %d, no. terms = %d  \n', size(M{1},1),size(N{1},1),length(M))
fprintf('\n No. equations (Ax=b): %d  \n', size(M{1},1)*size(N{1},1))

n=size(M{1},1); m=size(N{1},1);

% ----------------Choose shifts ------------------------------------
opts.tol=1e-4;  % NB: 'eigs' is VERY sensitive to this.

% This could be done more cheaply with data from a COARSER problem
% min and max eigs of hat{K}_{1} - 1st term could be strictly positive
% emin2 = 1e-6; 
emin2=1e-16; % eigs(M{2},M{1},1,'SA',opts); - for n>800 eigs doesnt work
emax2=eigs(M{2},M{1},1,'LA',opts);
n_m=size(M,2);
% emean=mean([(emin2),abs(emax2)]);
% alphas(1)=1-emean;
% alphas(2:n_m)=(1)*ones(n_m-1,1);        % center all other spectra around 1
    
%shifts = alphas'                       % Uncomment to print to screen

% Modify the orginal Ms and Ns (apply shifts) 
%  for k=2:size(M,2)
%      M{k}=M{k}+alphas(k-1)*M{1};
%      N{1}=N{1}-alphas(k-1)*N{k};
%  end
%  
% --------------Choose parameters s for MultiRB -------------------

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
% 
aa=emin2;%+alphas(1); 
bb=emax2;%+alphas(1);

% ZOLOTAREV POSITIVE IMAGINARY PARTS POLES
k = 4;      % rational degree
b = bb;     % sign function on [-10,-1]\cup [1,10]
r = rkfun.gallery('sign', k, b);
po = imag(poles(r));
s_parameter = po(po>=0 );

% OTHER POLES
% s_parameter = sqrt(aa*bb);
% s_parameter = complex(zeros(1, 6), logspace(log10(aa), log10(bb), 6));
% s_parameter = logspace(log10(aa),log10(bb),4)';

% GET_NODES2 POLES (SABINO THESIS)
% S_interval=[aa,bb];
% s_nodes = 4;                           % Choose 2 nodes (could vary)
% snew = get_nodes2(aa,bb,s_nodes);      % Use interval for A_1;
% s_parameter=snew;



% ------------ Options for MultiRB --------------------------------
res_method='4';                  % can't change this
rat_solve='1'; mmax=200;         % can change
param.max_space_dim=mmax;
param.period=1;                  % can change
param.rat_solve=rat_solve;
param.res_method=res_method;
tic;
fprintf('\n -------------------- MultiRB Solve ------------------------\n \n')

% [X1,X2,dimV,final_err,avg_inner,error_vec,iv_vec]=MultiRB_Poisson_rank1rhs_2sided(M,N,rhs1,rhs2,P,P1,param,s_parameter);
[X1,X2,dimV,final_err,avg_inner,error_vec,iv_vec]=MultiRB_noprec_Poisson_rank1rhs_2sided(M,N,rhs1,rhs2,param,s_parameter);
%X1(psort,:)=X1;
etoc=toc; 
fprintf('\n Total execution time: %9.4e seconds \n',etoc)
fprintf('\n ----------------------------------------------------------\n \n')
fprintf('no_terms   dim_M   dim_N   n_k   Rank    final_err   avg_inner  time(s)  \n')
fprintf('\n  %2d       %3d   %d    %3d    %2d    %9.4e    %4.2f     %9.4e  \n \n', [2, n, m, dimV, size(X1,2), final_err, avg_inner, etoc])
fprintf('\n----------------------------------------------------------\n \n')

% plot residual against iterations on semilogy plot
it = linspace(1, dimV, dimV);
semilogy(it, error_vec, 'x');
xlabel('Iterations');
ylabel('Log of the residual');hold on;
% upperbound_beckermann;
% plot(upper_bound, 'x'); hold off;