function [X1,X2,dimV,final_err,avg_inner,error_vec,iv_vec]=MultiRB_nonsymrhs(M,N,rhs1,rhs2,P,P1,param,snew)
%
% ---------------------- MultiRB solver -----------------------------------
%
% function [X1,X2,dimV,final_err,avg_inner,error_vec,iv_vec]=MultiRB(M,N,rhs1,rhs2,M,M1,param,snew);
%
% Solve  A x = b  arising from convection-diffusion using 
% matrix equation formulation:
%
%          M0 X N0 + M1 X N1 + .... + Mm X Nm = rhs1*rhs2'
%
% Note - problem is pre-processed (preconditioned + shifted) first by
% running SGFEM_matriver.m
%
% INPUTS:
% M            set of m+1 left-hand side matrices
% N            set of m+1 right-hand side matrices
% rhs1,rhs2    such that  b=vec(rhs1*rhs2')
% P, P1        preconditioner, P = P1 * P1'; set P1 to be identity
% param        contains MultiRB options, e,g:
%              param.max_space_dim=200;  (max approxn space dimension)
%              param.period=1;  (how often to check residual (no. its))
%              param.res_method='4'; (stopping criterion:  4 = relative diff soln)
%              param.rat_solve='1';  (RKS solver: 0=direct,1=iterative)
%
% OUTPUTS:
% X1, X2        X = X1*X2' final approximate solution
% dimV          dimension of final approximation space
% final_err     final relative difference in solution (see param.res_method)
% avg_inner     average no. inner iterations to solve projected problem

compute_period=param.period;
rat_solve=param.rat_solve;
res_method=param.res_method;
mmax=param.max_space_dim;
nterm=size(N,2);
tol=1e-5;       % Outer stopping tolerance (change if desired) 

tol_drop=.99;  % controls how many basis vectors to add at each iteration
nofirst=0;
Y=0;
Prhs1 = P1\rhs1; 
Prhs2 = P1\rhs2; 
normPrhs1=norm(Prhs1); normPrhs2 = norm(Prhs2);
% Include rhs in the basis.
V = (Prhs1)/normPrhs1;   
W = (Prhs2)/normPrhs2; 

tot_inner=[];
rhs1m=V'*Prhs1;
rhs2m=W'*Prhs2;
normR0=norm(full(Prhs1))*norm(full(Prhs2));
normR0_noprec=norm(full(rhs1))*norm(full(rhs2));
nrmres=1;
nrmres_noprec=1;
error_vec(1)=nrmres_noprec;

Y0=[];

s_nodes=length(snew);  % number of 's' parameters. 

% rat_space=1;
I=speye(size(M{1},1));
wrk1 = P1'\V;
wrk2 = P1'\W;
% here wrk1=wrk2

% Projected matrices -- M
for ind=1:nterm
    Mm{ind}=wrk2'*(M{ind}*wrk1); 
end

% Projected matrices -- N
for ind=1:nterm
    Nm{ind}=wrk2'*(N{ind}*wrk1);
end

m=size(N{1},1); n=size(M{1},1);
i=0; tot_it=0; ir=0; Y=[];

fprintf('      no.its      n_k       N_dim    error (rel diff)   no.inner its\n')

% main iteration
while (i < mmax & nrmres_noprec>tol) 
    i = i+1;
    tot_it = tot_it+1;
    if (i > size(V,2)), fprintf('exausted approx space\n'),end
    % Get new basis vectors - direct solver for solving all systems simultaneously
         ir = ir+1; if (ir>s_nodes),ir=1;end
            wrk1 = P1*V(1:n,i);
            for kk = 2:nterm
               v1(1:n,kk) = (M{kk}+snew(ir)*M{1})\wrk1;
            end
            v1 = P1'*v1;
            
            wrk2 = P1*W(1:m,i);
            for kk = 2:nterm
               v2(1:m,kk) = (N{kk}+snew(ir)*N{1})\wrk2;
            end
            v2 = P1'*v2;

    v1 = v1 - V*(V'*v1); v1 = v1 - V*(V'*v1);
    v2 = v2 - W*(W'*v2); v2 = v2 - W*(W'*v2);
    
    % Deselect new basis vectors (TRUNCATE)
    [uu1,ss1,vv1]=svd(v1,0);
    ss1=diag(ss1);
    
    [uu2,ss2,vv2]=svd(v2,0);
    ss2=diag(ss2);
    
    iv=size(V,2);
    iv_vec(i)=iv;
   
    if ss1(1,1)>1e-12 % ss1(1,1)  = 7.9603 e-17 for n=5
        addv=1; 
        l1=cumsum(ss1)/sum(ss1); il1=find(l1>=tol_drop,1); 
        vnew1=uu1(:,1:il1);
        Pvnew1=P1'\vnew1;
        V(1:n,iv+1:iv+il1)=vnew1; % increase the space V        
    else
       addv=0;
    end
    
    iw=size(W,2);
    iw_vec(i)=iw;
    if ss2(1,1)>1e-12
        addv=1; 
        l2=cumsum(ss2)/sum(ss2); il2=find(l2>=tol_drop,1); 
        vnew2=uu2(:,1:il2);
        Pvnew2=P1'\vnew2;
        W(1:m,iw+1:iw+il2)=vnew2; % increase the space W       
    else
       addv=0;
    end

    iv=size(V,2);
    iw=size(W,2);
  
    vector_length=iv*iw;
    
    % Expand projected matrices -- M
    Mm{1}=speye(iv);
    ivnew=size(vnew1,2);

    for ind=2:nterm
        wrk1=M{ind}*Pvnew1;
        Pwrk1=P1\wrk1;
        newk1=V(1:n,1:(iv-ivnew))'*Pwrk1; 
        Mm{ind}(1:iv-ivnew,iv-ivnew+1:iv)=newk1;  
        Mm{ind}(iv-ivnew+1:iv,1:iv-ivnew)=newk1';
        Mm{ind}(iv-ivnew+1:iv,iv-ivnew+1:iv)=Pwrk1'*V(:,iv-ivnew+1:iv);
    end

    if addv, rhs1m=[rhs1m; vnew1'*Prhs1];end

    % Expand projected matrices -- N
%     Nm{1}=speye(iv);
    iwnew=size(vnew2,2);

    for ind=1:nterm
        wrk2=N{ind}*Pvnew2;
        Pwrk2=P1\wrk2;
        newk2=V(1:m,1:(iw-iwnew))'*Pwrk2; 
        Nm{ind}(1:iw-iwnew,iw-iwnew+1:iw)=newk2;  
        Nm{ind}(iw-iwnew+1:iw,1:iw-iwnew)=newk2';
        Nm{ind}(iw-iwnew+1:iw,iw-iwnew+1:iw)=Pwrk2'*W(:,iw-iwnew+1:iw);
    end

    if addv, rhs2m=[rhs2m; vnew2'*Prhs2];end

    % Periodically compute the approx solution and residual
    if (rem(tot_it,compute_period)==0)
%         Yk = (kron(Nm{1}', Mm{1}) + kron(Nm{2}', Mm{2}) + kron(Nm{3}', Mm{3}) + kron(Nm{4}', Mm{4}))\(kron(rhs2m', rhs1m));
        tol_inner=nrmres_noprec*1e-5;
        y0=zeros(iv,iw);
        if (nofirst) 
            y0(1:size(Y,1),1:size(Y,2))=Y;
        end

        [Y,iteraY]=cgkron(Mm,Nm,rhs1m*rhs2m',y0,iv*iw,tol_inner);  %,iv,iw); % inner solver: cg

        tot_inner=[tot_inner,iteraY];
        nofirst=1;
        % Compute relative variation in the solution
        Yer=Y; [n0,m0]=size(Y0); Yer(1:n0,1:m0)=Yer(1:n0,1:m0)-Y0; 
        ErrY=norm(Yer,'fro'); %/normR0;
        nrmres_noprec=ErrY/norm(Y,'fro'); nrmres=nrmres_noprec;
        Y0=Y;
        
        % Print progress to screen
         fprintf('\n        %2d        %3d       %d      %8.6e         %2d \n', [i,iv,iw,nrmres,iteraY])
    end

   error_vec(i+1)=nrmres_noprec;
end

fprintf('\n Total iterations: %d \n\n', i)

%------------------------------------------------
%
% Optional - Estimate Rank + Compute Factors
% 
%-----------------------------------------------

% Y = reshape(Y, n, n);
[uu,ss,vv]=svd(Y,0);
ns=sum(diag(ss)/ss(1,1)>tol/iv);

% Determine the solution factors
if (size(V,2)<size(W,1))
   X1=P1'\(V*uu(:,1:ns)); 
   X2=(ss(1:ns,1:ns)*vv(:,1:ns)'*W')';
else
   X1=P1'\(V*(uu(:,1:ns)*ss(1:ns,1:ns))); 
   X2=W*vv(:,1:ns);
end

dimV=size(V,2); final_err=nrmres; avg_inner=mean(tot_inner);

