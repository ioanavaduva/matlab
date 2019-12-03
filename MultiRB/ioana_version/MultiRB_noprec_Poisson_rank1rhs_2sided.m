function [X1,X2,dimV,final_err,avg_inner,error_vec,iv_vec]=MultiRB_noprec_Poisson_rank1rhs_2sided(M,N,rhs1,rhs2,param,snew)
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
tol=1e-9;       % Outer stopping tolerance (change if desired) 

tol_drop=.99;  % controls how many basis vectors to add at each iteration
nofirst=0;
Y=0;
normrhs1=norm(rhs1); normrhs2 = norm(rhs2);
% Include rhs in the basis.
V = (rhs1)/normrhs1;   
W = (rhs2)/normrhs2; 
% at this stage V and W are the same

tot_inner=[];
rhs1m=V'*rhs1;
rhs2m=W'*rhs2;
% here rhs1m=rhs2m
% normR0=norm(full(rhs1))*norm(full(rhs2));
% normR0_noprec=norm(full(rhs1))*norm(full(rhs2));
% nrmres=1;
nrmres_noprec=1;
error_vec(1)=nrmres_noprec;

Y0=[];

s_nodes=length(snew);  % number of 's' parameters. 

% rat_space=1;
% I=speye(size(M{1},1));
wrk1 = V;
wrk2 = W;
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
            wrk1 = V(1:n,i);
            for kk = 2:nterm
               v1(1:n,kk) = (M{kk}+snew(ir)*M{1})\wrk1;
            end
%             size(m)
%             size(n)
            wrk2 = wrk1; 
% % %             wrk2 = P1*W(1:m,i)
%             for kk = 2:nterm
%                v2(1:m,kk) = (M{kk}+snew(ir)*M{1})\wrk2;
%             end

    v1 = v1 - V*(V'*v1); v1 = v1 - V*(V'*v1);
v2=v1;
%     v2 = v2 - W*(W'*v2); v2 = v2 - W*(W'*v2)
    
    % Deselect new basis vectors (TRUNCATE)
    [uu1,ss1,vv1]=svd(v1,0);
    ss1=diag(ss1);
    
    [uu2,ss2,vv2]=svd(v2,0);
    ss2=diag(ss2);
    
% separate computations
%     iv=size(V,2); 
%     iv_vec(i)=iv;
%  
%     if ss1(1,1)>1e-12
%         addv=1; 
%         l1=cumsum(ss1)/sum(ss1); il1=find(l1>=tol_drop,1); 
%         vnew1=uu1(:,1:il1);
%         Pvnew1=P1'\vnew1;
%         V(1:n,iv+1:iv+il1)=vnew1; % increase the space V        
%     else
%        addv=0;
%     end
%     
%     iw=size(W,2); 
%     iw_vec(i)=iw;
%     
%     if ss2(1,1)>1e-12
%         addv=1; 
%         l2=cumsum(ss2)/sum(ss2); il2=find(l2>=tol_drop,1); 
%         vnew2=uu2(:,1:il2);
%         Pvnew2=P1'\vnew2;
%         W(1:m,iw+1:iw+il2)=vnew2; % increase the space W       
%     else
%        addv=0;
%     end

% combined if ss >1e-12
    iv=size(V,2); 
    iv_vec(i)=iv;
    iw=size(W,2); 
    iw_vec(i)=iw;
%     sone=ss1(1,1);
%     stwo=ss2(1,1);
    if ss1(1,1)>1e-12 && ss2(1,1)>1e-12
        addv=1; 
        l1 = cumsum(ss1)/sum(ss1); il1=find(l1>=tol_drop,1); 
        vnew1 = uu1(:,1:il1);
        Pvnew1 = vnew1;
        V(1:n,iv+1:iv+il1) = vnew1;
        l2 = cumsum(ss2)/sum(ss2); il2 = find(l2>=tol_drop,1); 
        vnew2 = uu2(:,1:il2);
        Pvnew2 = vnew2;
%         W(1:m,iw+1:iw+il2)=vnew2;
        W=V;
    else
        addv=0;
    end
    
%  V
%  W
% W'*V
% V'*V
    iv=size(V,2);
    iw=size(W,2);
  
    vector_length=iv*iw;
    
    % Expand projected matrices -- M
    Mm{1}=speye(iv);
    ivnew=size(vnew1,2);


    for ind=2:nterm
        wrk1=M{ind}*Pvnew1;
        Pwrk1 = wrk1;
        newk1 = V(1:n,1:(iv-ivnew))'*Pwrk1; 
        Mm{ind}(1:iv-ivnew,iv-ivnew+1:iv) = newk1;  
        Mm{ind}(iv-ivnew+1:iv,1:iv-ivnew) = newk1';
        Mm{ind}(iv-ivnew+1:iv,iv-ivnew+1:iv) = Pwrk1'*V(:,iv-ivnew+1:iv);
    end

    if addv, rhs1m=[rhs1m; vnew1'*rhs1];end

    % Expand projected matrices -- N
%     Nm{1}=speye(iv);
    iwnew = size(vnew2,2);


    for ind = 1:nterm
        wrk2 = N{ind}*Pvnew2;
        Pwrk2 = wrk2;
        newk2 = V(1:m,1:(iw-iwnew))'*Pwrk2; 
        Nm{ind}(1:iw-iwnew,iw-iwnew+1:iw) = newk2;  
        Nm{ind}(iw-iwnew+1:iw,1:iw-iwnew) = newk2';
        Nm{ind}(iw-iwnew+1:iw,iw-iwnew+1:iw) = Pwrk2'*W(:,iw-iwnew+1:iw);
    end % 28/10/19 Nm{2} should stay identity? for 5x5 pb now fist entry is -0.8944
    
    if addv, rhs2m = [rhs2m; vnew2'*rhs2];end
%     size(rhs1m*rhs2m')
%     size(kron(Nm{1}', Mm{1}))

    % Periodically compute the approx solution and residual
    if (rem(tot_it,compute_period)==0)
%         rhsm = rhs1m*rhs2m';
%         rhsmm = rhsm(:);
%         Y = lyap(Nm{1}, rhsm);
%         Y = (kron(Nm{1}', Mm{1}) + kron(Nm{2}', Mm{2}))\rhsmm; %+ kron(Nm{3}', Mm{3}) + kron(Nm{4}', Mm{4}))\(kron(rhs2m', rhs1m));
        tol_inner = nrmres_noprec*1e-1;
        y0 = zeros(iv,iw);
        if (nofirst) 
            y0(1:size(Y,1),1:size(Y,2)) = Y;
        end
%         size(Mm{1})
%         size(Mm{2})
%         size(rhs1m*rhs2m')
       [Y,iteraY] = cgkron(Mm,Nm,rhs1m*rhs2m',y0,iv*iw,tol_inner); % inner solver: cg

%  Y = lyap(Mm{2}, rhs1m*rhs2m');
%         iteraY = 1;
        tot_inner = [tot_inner,iteraY];
%         nofirst = 1;
        
% Compute the residual - exactly
        Yer = Y; [n0,m0] = size(Y0); Yer(1:n0,1:m0) = Yer(1:n0,1:m0)-Y0; 
% size(V)
% size(Y)
% size(W')
        X_hat = V*Y*W';
% %         ErrY=norm(Yer,'fro'); %/normR0;
        ErrY = norm(X_hat*M{2} + M{2}*X_hat - rhs1*rhs2', 'fro');
% %         nrmres_noprec=ErrY/norm(Y,'fro'); nrmres=nrmres_noprec;
        nrm_rhs = norm(rhs1*rhs2', 'fro');
        nrmres_noprec = ErrY/nrm_rhs;
%         Y0 = Y;

% Compute the residual - as in Druskin, Simoncini
% % computed residual   (exact, in exact arithmetic)
%      u1=newAv-VV(1:n,1:js)*g;
%      d=-VV(1:n,1:js)*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
%      U=[-V*s(end),  d u1 ];
%      rr=qr(full(U),0); rr=triu(rr(1:size(rr,2),:));
%      nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/(nrmb+singE*nrma*nrmx);
%      nrmrestot=[nrmrestot,nrmres];
% 
%      disp([i,nrmres]) 
%      if (nrmres<tol), break,end

% Compute the residual - as in Kirsten, Simoncini

% % computed residual   (exact, in exact arithmetic) cheaper computation possible
%      u1=newAv-VV*g;
%      d=-VV*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
%      U=full([-V*s(end),  d u1 ]);
%      rr=qr(U,0); rr=triu(rr(1:size(rr,2),:));
%      
% 
% %backward error
%      nrmres = norm(rr*sparse([O I O; I O I; O I O ])*rr','fro');
% 
% % relative residual norm
%      nrmresnew = (nrmres)/nrmb;
%  
%      nrmrestotnew = [nrmrestotnew, nrmresnew];
% 
%      dim = size(VV,2);
%      dim1 = [dim1,dim];
% 
%      disp([i,nrmresnew])
% 
%      if (nrmresnew<tol), 
%         break
%      end
%   end 
        
        % Print progress to screen
         fprintf('\n        %2d        %3d       %d      %8.6e         %2d \n', [i,iv,iw,nrmres_noprec,iteraY])
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
[uu,ss,vv] = svd(Y,0);
ns = sum(diag(ss)/ss(1,1)>tol/iv);

% Determine the solution factors
if (size(V,2)<size(W,1))
   X1 = V*uu(:,1:ns); 
   X2=(ss(1:ns,1:ns)*vv(:,1:ns)'*W')';
else
   X1 = V*(uu(:,1:ns)*ss(1:ns,1:ns)); 
   X2=W*vv(:,1:ns);
end

dimV = size(V,2); final_err = nrmres_noprec; avg_inner = mean(tot_inner);

