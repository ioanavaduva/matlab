% Inhomogeneous Helmhotz equation: we are solving the generalized Sylvester equation 
%
% A*X+X*B^T+ N*X*N^T= C*C^T (**)
%
% by the extended block Krylov subspace method (EKSM) where the initial
% block is selected exploiting the low-rank commutation property of the
% coefficients.
%
% Equation (**) comes from the discretization by finite differences of
% the inhomogeneous Helmhotz equation 
%   -u_xx-u_yy+chi(x,y)u=f(x,y)
%   u(x,0)=0;   u(x,1)=0;   homogeneous Dirichlet boundary conditions
%   u(x,y+1)=u(x,y);        periodic boundary conditions 
%
%   chi(x,y) is the indicator function of the region (x>1/2)(y>1/2) 
%   f(x,y) is a source term, localized in a square
%    
% For more details see Example 4.3 in:
%
% Elias Jarlebring, Giampaolo Mele, Davide Palitta and Emil Ringh
% Krylov methods for low-rank commuting generalized Sylvester equations
% April 2017, ArXiv: 1704.02167
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

clear all

n_vect=[5e3 1e4 5e4 1e5 5e5 1e6];
n_vect=100
for i=1:length(n_vect)    
    n=n_vect(i);

    % Generate data    
    C=10*ones(n,1);
    C(1:n/4)=0; C(n/2:end)=0;

    e = ones(n,1);
    B = -spdiags([e -2*e e], -1:1, n, n);
    A = -spdiags([e -2*e e], -1:1, n, n);
    A(1,end)=-1;    A(end,1)=-1;
    B=n^2*B;        A=n^2*A;
        
    % compute the other coefficients of the matrix equation
    O=sparse(n/2,n/2);
    I=speye(n/2,n/2);
    N=[O, O; O, I];
    
    % low-rank commutating factors. 
    I=speye(n);
    v=n*I(:,n/2+1);
    w=n*I(:,n/2);
    e1=n*I(:,1);  
    eN=n*I(:,end);
    % REMARK: It holds 
    % A*N-N*A = [v,w,e1,en]*[w,-v,-en,e1]^T
    % B*N-N*B = [v,w]*[w,-v]^T
    norm(A*N-N*A-[v,w,e1,eN]*[w,-v,-eN,e1]',1)
    norm(B*N-N*B-[v,w]*[w,-v]',1)
   
   pause
    % A is singular, so we shift it AA=A+I.
    % we thus solve the equivalent problem: 
    % (A+I)*X + X*B^T + N*X*N^T -X = C*C^T
    % the space we need to generate is the same, the projected problem will
    % be different.
    N_cell=cell(2);
    N_cell{1}=N;
    N_cell{2}=I;
    M_cell=N_cell;
    M_cell{2}=-I;
    
    U=full([v,w,e1,eN]);
    Q=full([v,w]);
    
    starting_block_left=[C,U];
    starting_block_right=[C,Q];
    
    tol=-Inf;
    % Precompute the LU factors
    AA=A+I;
    [LA,UA]=lu(AA);
    [LB,UB]=lu(B);
    tolY=1e-12;
    maxit=30;
    

    [X1,X2,res]=EKSM_genSylv_LRcomm(AA,LA,UA,B,LB,UB,N_cell,starting_block_left,M_cell,starting_block_right,C,C,maxit,tol,tolY);

    % for not too large problems you can check the relative residual norm 
    % norm((A*X1)*X2'+X1*(X2'*B')+(N*X1)*(X2'*N')-C*C','fro')/norm(C*C','fro')
    pause

end

