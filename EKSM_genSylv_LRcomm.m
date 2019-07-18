function [X1,X2,res]=EKSM_genSylv_LRcomm(A,LA,UA,B,LLB,UB,N,starting_block_left,M,starting_block_right,C1,C2,maxit,tol,tolY)
% function [X1,X2,res]=EKSM_genSylv_LRcomm(A,LA,UA,B,LB,UB,N,starting_block_left,M,starting_block_right,C1,C2,maxit,tol,tolY)
% Extended Krylov subspace method for generalized Sylvester equation with
% low-rank commuting coefficients.
%
% We solve
%
% A X + X B' + \sum_{i=1}^m N_i*X*M_i^T = C_1*C_2^T
%
% assuming the Sylvester part being the dominant one and 
%
% [A,N_i]=A*N_i-N_i*A=U*\tilde U^T for all i=1,...,m,
% [B,M_i]=B*M_i-M_i*B=Q*\tilde Q^T for all i=1,...,m,
%
% where U, \tilde U, Q, \tilde Q are all tall matrices.
%
% Input:
% A = LA*UA;   LU factorization of A, 
% B = LB*UB;   LU factorization of B, 
% N: cell that contains the m matrices N_i, N{i}=N_i
% starting_block_left: low-rank matrix we use to build the left space:
%                 EK_k^\square(A,starting_block_left)
% M: cell that contains the m matrices M_i, M{i}=M_i
% starting_block_right: low-rank matrix we use to build the right space:
%                 EK_k^\square(B,starting_block_right)
% C1,C2: low-rank factors of the right-hand side
% maxit = max space dimension
% tol = max final accuracy (in terms of relative residual)
% tolY = at convergence a truncated SVD up to tolY is performed in order to further reduce the rank of the 
%         computed solution 
%
% Output:
% X1,X2: solution factors   X \approx X1*X2^T
% res:  history of residual norms
%
% REMARK: the Matlab Control Toolbox, and in particular the function
% lyap.m, is requested.
%
% Please contact D. Palitta for any problem you may encouter when running the code.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
% FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
% IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

tic;

% check the input
if iscell(A) || iscell(C1) || iscell(B) || iscell(C2)
    error('A, B, C1 and C2 must all be defined as matrices')
elseif ~iscell(N) || ~iscell(M)
    error('N and M must all be defined as cells')
end


nrmC1=norm(C1,'fro');  
nrmC2=norm(C2,'fro');  
normC=nrmC1*nrmC2;

% size of the "real" rhs
ss1=size(C1,2);
ss2=size(C2,2);

m=length(N);


[n,sh]=size(starting_block_left);
[nB,sh1]=size(starting_block_right);


odds=[];
oddsB=[];

% if A and B are symmetric the Lanczos method can be employed in the basis
% construction setting kA_max=kB_max=3, but the Arnoldi procedure is much more
% stable
kA_max=maxit;
kB_max=maxit;
sA=2*sh;
sB=2*sh1;

% Sequence in A
starting_block_leftA=UA\(LA\starting_block_left);  
[U(1:n,1:sA),beta]=qr([starting_block_left,starting_block_leftA],0);
gamma1=beta(1:ss1,1:ss1);
TTN=cell(m);
TTN_old=cell(m);
ibeta=inv(beta);
H=zeros((maxit+1)*sA,maxit*sA); L=zeros((maxit+1)*sA,maxit*sA);
T=zeros(maxit*sA,maxit*sA);

% Sequence in B'
starting_block_rightB=UB\(LLB\starting_block_right);
[W(1:nB,1:sB),betaB]=qr([starting_block_right,starting_block_rightB],0);
theta1=betaB(1:ss2,1:ss2);
TTM=cell(m);
TTM_old=cell(m);
ibetaB=inv(betaB);
HB=zeros((maxit+1)*sB,maxit*sB); LB=zeros((maxit+1)*sB,maxit*sB);
TB=zeros(maxit*sB,maxit*sB);


beta2=gamma1*theta1';
t1=toc;
time_orth=0;
time_inner=0;
time_rest=0;
res=zeros(maxit+2,1);

for j=1:maxit,

    tic
    jms=(j-1)*sA+1;j1s=(j+1)*sA;js=j*sA;js1=js+1; jsh=(j-1)*sA+sh;

    % Sequence in A
    Up(1:n,1:sh)   = A*U(1:n,jms:jsh); 
    Up(1:n,sh+1:sA) = UA\(LA\U(1:n,jsh+1:js));

    %new bases block (modified gram)
    for l=1:2
        k_min=max(1,j-kA_max);
        for kk=k_min:j
            k1=(kk-1)*sA+1; k2=kk*sA;
            coef= U(1:n,k1:k2)'*Up;
            H(k1:k2,jms:js) = H(k1:k2,jms:js)+ coef; 
            Up = Up - U(:,k1:k2)*coef; 
        end
    end
    
    % normalization
    if (j<=maxit)
       [U(1:n,js1:j1s),H(js1:j1s,jms:js)]=qr(Up,0);
       hinv=inv(H(js1:j1s,jms:js));
    end
    
    t2=toc;
    time_orth=time_orth+t2;
    
    tic
    
    I=speye(js+sA);
    if (j==1),
      L(1:j*sA+sh,(j-1)*sh+1:j*sh) =...
      [ H(1:sA+sh,1:sh)/ibeta(1:sh,1:sh), speye(sA+sh,sh)/ibeta(1:sh,1:sh)]*ibeta(1:sA,sh+1:sA);
    else
      L(1:j*sA+sA,(j-1)*sh+1:j*sh) = L(1:j*sA+sA,(j-1)*sh+1:j*sh) + H(1:j*sA+sA,jms:jms-1+sh)*rho;
    end
    odds  = [odds, jms:(jms-1+sh)];   % store the odd block columns
    evens = 1:js; evens(odds)=[];
    T(1:js+sA,odds)=H(1:js+sA,odds);   %odd columns

    T(1:js+sh,evens)=L(1:js+sh,1:j*sh);   %even columns
    L(1:j*sA+sA,j*sh+1:(j+1)*sh) = ...
       ( I(1:j*sA+sA,(js-sh+1):js)- T(1:js+sA,1:js)*H(1:js,js-sh+1:js))*hinv(sh+1:sA,sh+1:sA);
    rho = hinv(1:sh,1:sh)\hinv(1:sh,sh+1:sA);
    
    jmsA=jms;j1sA=j1s;jsA=js;js1A=js1; jshA=jsh;

    tA=toc;
    
    %Sequence in B'
    tic
    jms=(j-1)*sB+1;j1s=(j+1)*sB;js=j*sB;js1=js+1; jsh=(j-1)*sB+sh1;

    Wp(1:nB,1:sh1)   = B*W(1:nB,jms:jsh); 
    Wp(1:nB,sh1+1:sB) = UB\(LLB\W(1:nB,jsh+1:js));

    %new bases block (modified gram)
    for l=1:2
        k_min=max(1,j-kB_max);
        for kk=k_min:j
            k1=(kk-1)*sB+1; k2=kk*sB;
            coef= W(1:nB,k1:k2)'*Wp;
            HB(k1:k2,jms:js) = HB(k1:k2,jms:js)+ coef; 
            Wp = Wp - W(:,k1:k2)*coef; 
        end
    end
    
    % normalization 
    if (j<=maxit)
       [W(1:nB,js1:j1s),HB(js1:j1s,jms:js)]=qr(Wp,0);       
       hinvB=inv(HB(js1:j1s,jms:js));
    end
    
    t2=toc;
    time_orth=time_orth+t2;
    
    tic
    
    I=speye(js+sB);
    if (j==1),
      LB(1:j*sB+sh1,(j-1)*sh1+1:j*sh1) =...
      [ HB(1:sB+sh1,1:sh1)/ibetaB(1:sh1,1:sh1), speye(sB+sh1,sh1)/ibetaB(1:sh1,1:sh1)]*ibetaB(1:sB,sh1+1:sB);
    else
      LB(1:j*sB+sB,(j-1)*sh1+1:j*sh1) = LB(1:j*sB+sB,(j-1)*sh1+1:j*sh1) + HB(1:j*sB+sB,jms:jms-1+sh1)*rhoB;
    end
    oddsB = [oddsB, jms:(jms-1+sh1)];   % store the odd block columns
    evens = 1:js; evens(oddsB)=[];
    TB(1:js+sB,oddsB)=HB(1:js+sB,oddsB);   %odd columns

    TB(1:js+sh1,evens)=LB(1:js+sh1,1:j*sh1);   %even columns
    LB(1:j*sB+sB,j*sh1+1:(j+1)*sh1) = ...
       ( I(1:j*sB+sB,(js-sh1+1):js)- TB(1:js+sB,1:js)*HB(1:js,js-sh1+1:js))*hinvB(sh1+1:sB,sh1+1:sB);
    rhoB = hinvB(1:sh1,1:sh1)\hinvB(1:sh1,sh1+1:sB);

    k=j;

    tB=toc;
    
    tic
    
    % Solve reduced multiterm eqn avoiding explicit projections
    for i=1:m
        % compute the projection of N{i} onto new basis blocks
        vec1=U(1:n,1:jsA)'*(N{i}*U(1:n,jmsA:jsA));
        vec2=(U(1:n,jmsA:jsA)'*N{i})*U(1:n,1:jmsA-1);
        % assemble with the old projected N_i
        TTN_temp=zeros(jsA);
        TTN_temp(1:jmsA-1,1:jmsA-1)=TTN_old{i};
        TTN_temp(:,jmsA:jsA)=vec1;
        TTN_temp(jmsA:jsA,1:jmsA-1)=vec2;
        TTN{i}=TTN_temp;
        TTN_old{i}=TTN{i};
        
        % compute the projection of M{i} onto the new basis blocks
        vec1=W(1:nB,1:js)'*(M{i}*W(1:nB,jms:js));
        vec2=(W(1:nB,jms:js)'*M{i})*W(1:nB,1:jms-1);
        % assemble with the old projected E_i
        TTM_temp=zeros(js);
        TTM_temp(1:jms-1,1:jms-1)=TTM_old{i};
        TTM_temp(:,jms:js)=vec1;
        TTM_temp(jms:js,1:jms-1)=vec2;
        TTM{i}=TTM_temp;
        TTM_old{i}=TTM{i};
        
    end
    
    t_proj=toc;
   
    tic
    Y=GenSylv_smallscale(T(1:jsA,1:jsA),TB(1:js,1:js)',TTN,TTM,eye(k*sA,ss1)*beta2*eye(k*sB,ss2)',1e-16,100);
    t4=toc;
    time_inner=time_inner+t4;
    
    tic
    
    cc  = [H(js1A:j1sA,jsA-sA+1:jsA-sh), L(js1A:j1sA,(j-1)*sh+1:j*sh)];
    ccB = [HB(js1:j1s,js-sB+1:js-sh), LB(js1:j1s,(j-1)*sh+1:j*sh)];
    
    % compute the residual norm
    res(k)=sqrt(norm(cc*Y(jsA-sA+1:jsA,:),'fro')^2+norm(Y(:,js-sB+1:js)*ccB','fro')^2)/normC;
    t_res=toc;
    
    time_rest=time_rest+tA+tB+t_proj+t_res;
    
if (res(k)<tol), break, end
end

tic
% Reduce rank of solution, if needed
[uY,sY,vY] = svd(Y); 
[sY,id] = sort(diag(sY));
sY = flipud(sY);
uY = uY(:,id(end:-1:1));
vY = vY(:,id(end:-1:1));
is = sum(abs(sY)>tolY);
Y1 = uY(:,1:is)*diag(sqrt(sY(1:is))); 
Y2 = vY(:,1:is)*diag(sqrt(sY(1:is)));

X1=U(1:n,1:jsA)*Y1; X2=W(1:nB,1:js)*Y2;

t6=toc;
time_rest=time_rest+t6+t1;


time_tot=time_orth+time_inner+time_rest;
time_orth_perc=time_orth*100/time_tot;
time_inner_perc=time_inner*100/time_tot;
time_rest_perc=time_rest*100/time_tot;

fprintf('Its: %d Left Space dim %d Right Space dim %d Solution rank %d \n',k,jsA,js,is);
fprintf('Time_orth: %10.5e Time_inner: %10.5e Time_rest: %10.5e Time_tot: %10.5e \n',time_orth,time_inner,time_rest,time_tot)
fprintf('Time_orth_perc: %10.5e Time_inner_perc: %10.5e Time_rest_perc: %10.5e \n',time_orth_perc,time_inner_perc,time_rest_perc)

res=res(1:k);

return
