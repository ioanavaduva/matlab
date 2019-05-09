function [Z,er2]=kpik(A,E,LE,B,m,tol,tolY)
%function [Z,er2]=kpik(A,E,LE,B,m,tol,tolY)
%
% Approximately solve
%
%       A X E + E X A' + BB' = 0
%
% by means of the extended Krylov subspace method
%
% Author of the matlab code:  Valeria Simoncini
% version 1.0
%  The code uses the function lyap.m of the Control Matlab Toolbox
%
%
% Input
%  A   coeff matrix, A < 0
%  E   coeff matrix, spd, 
%  LE  lower triang factor of coeff matrix
%  B   factor of rhs,   nxk matrix with k << n
%  m   max space dimension, say sqrt(size(A))
%  tol stopping tolerance, with stopping criterion
%          ||LE\A X LE  + LE' X A'/LE'-LE\BB'/LE'||
%          ----------------------------------------  < tol
%      ||LE\BB'/LE'|| + ||E^{-1}|| ||A|| ||LE'X LE ||
%      computed in a cheap manner
%
%  Output:
%  Z   solution factor   X = Z Z'
%  er2 history of scaled residual, as above
%
%
% Comments: 
% * The projected solution is computed at each iteration
%   As an alternative, a periodic computation could be considered.
% * This code performs a factorization of A. As an alternative,
%   iterative solves could be considered.
% 
%  Please contact V. Simoncini for any problem you may encouter when
%  running the code
%
% If you use this code, please cite the following article:
%
% V. Simoncini
% A new iterative method for solving large-scale Lyapunov matrix equations, 
% SIAM J.  Scient. Computing, v.29, n.3 (2007), pp. 1268-1288.
%
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
%FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
%COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
%IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
%CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
%

tic;

rhs=LE\B;
nrmb=norm(rhs,'fro')^2;  
nrma=norm(A,'fro');
sqrt2=sqrt(2);
er2=zeros(m,1);

[n,sh]=size(rhs);

Y=[];
odds=[];

if (norm(E-speye(n),'fro')>1e-14)
  condestE=condest(E);
  singE=condestE/norm(E,'fro');
else
  singE=1;
end


if norm(A-A',1)<1e-14,
     UA = chol(-A); LA = -UA';
     fprintf('A sym. Completed Chol factorization\n')
     k_max =2;
else
     [LA,UA]=lu(A);
     fprintf('A nonsym. Completed LU factorization\n')
     k_max = m;
end

s=2*sh;
rhs1=LE'*(UA\(LA\(LE*rhs)));

[U(1:n,1:s),beta]=qr([rhs,rhs1],0);
ibeta=inv(beta(1:s,1:s));
beta = beta(1:sh,1:sh); 
beta2=beta*beta';    
H=zeros((m+1)*s,m*s);
T=zeros((m+1)*s,m*s);
L=zeros((m+1)*s,m*s);
fprintf('      it        backward err\n')

for j=1:m,

    jms=(j-1)*s+1;j1s=(j+1)*s;js=j*s;js1=js+1; jsh=(j-1)*s+sh;
    Up(1:n,1:sh) = LE\(A*(LE'\U(:,jms:jsh))); 
    Up(1:n,sh+1:s) = LE'*(UA\(LA\(LE*U(1:n,jsh+1:js))));

%new bases block (modified gram)
    for l=1:2
        k_min=max(1,j-k_max);
        for kk=k_min:j
            k1=(kk-1)*s+1; k2=kk*s;
            coef= U(1:n,k1:k2)'*Up;
            H(k1:k2,jms:js) = H(k1:k2,jms:js)+ coef;
            Up = Up - U(:,k1:k2)*coef; 
        end
    end

  if (j<=m)
    [Up,H(js1:j1s,jms:js)] = qr(Up,0);
    hinv=inv(H(js1:j1s,jms:js));
  end

  I=eye(js+s);
  if (j==1),
      L(1:j*s+sh,(j-1)*sh+1:j*sh) =...
      [ H(1:s+sh,1:sh)/ibeta(1:sh,1:sh), speye(s+sh,sh)/ibeta(1:sh,1:sh)]*ibeta(1:s,sh+1:s);
  else
      L(1:j*s+s,(j-1)*sh+1:j*sh) = L(1:j*s+s,(j-1)*sh+1:j*sh) + H(1:j*s+s,jms:jms-1+sh)*rho;
  end

  odds = [odds, jms:(jms-1+sh)];   % store the odd block columns
  evens = 1:js; evens(odds)=[];
  T(1:js+s,odds)=H(1:js+s,odds);   %odd columns

  T(1:js+sh,evens)=L(1:js+sh,1:j*sh);   %even columns
  L(1:j*s+s,j*sh+1:(j+1)*sh) = ...
       ( I(1:j*s+s,(js-sh+1):js)- T(1:js+s,1:js)*H(1:js,js-sh+1:js))*hinv(sh+1:s,sh+1:s);
  rho = hinv(1:sh,1:sh)\hinv(1:sh,sh+1:s);

  k=j;
  Y = lyap((T(1:js,1:js)),eye(k*s,sh)*beta2*eye(k*s,sh)');
% safeguard to preserve symmetry
  Y = (Y+Y')/2;

  cc = [H(js1:j1s,js-s+1:js-sh), L(js1:j1s,(j-1)*sh+1:j*sh)];
  nrmx=norm(Y,'fro');
  er2(k)=sqrt2*norm(cc*Y(js-s+1:js,:),'fro')/(nrmb+singE*nrma*nrmx);

  disp([k,er2(k)])

  if (er2(k)<tol), 
     break
  else
     U(1:n,js1:j1s)=Up;
  end
end

% Done
% reduce solution rank if needed
  [uY,sY]=eig(Y); [sY,id]=sort(diag(sY));
  sY=flipud(sY); uY=uY(:,id(end:-1:1));
  is=sum(abs(sY)>tolY);
  Y0 = uY(:,1:is)*diag(sqrt(sY(1:is)));
  Z = LE'\(U(1:n,1:js)*Y0);
  er2=er2(1:j);
  total_time=toc;
  fprintf('   its           comp.res.   space dim.   CPU Time\n')
  disp([k,er2(j),js,total_time])

return
