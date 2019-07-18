function [X,solves,res,err]=GenSylv_smallscale(A,B,N,M,C,tol,max_solves,Xexact)
% function [X,solves,res,err]=GenSylv_smallscale(A,B,N,M,C,tol,max_solves,Xexact)
% Function that approximately solves the small-scale matrix equation
%
% A*X + X*B^T + \sum_{i=1}^m N_i*X*M_i^T = C (**).
%
% The solution X to (**) is approximated as
%   X = \sum_{k=0}^\ell Y_k
% where
% Y_0 is such that  A * Y_0 + Y_0 * B^T = C
% Y_j, j>0, is such that A * Y_j + Y_j * B^T = -\sum_{i=1}^m N_i*Y_{j-1}*M_i^T
%
% Elias Jarlebring, Giampaolo Mele, Davide Palitta and Emil Ringh
% Krylov methods for low-rank commuting generalized Sylvester equations
% April 2017, ArXiv: 1704.02167
%
% INPUT:
% A, B, N, M: coefficient matrices. A and B must be passed as matrices
%             while N and M must be cells that collect the matrices N_i, 
%             M_i, that is N{i}=N_i and M{i}=M_i
% C: given term
% tol: tolerance
% max_solves: maximum number of Sylvester solves to perform
% Xexact: exact solution to (**), if available
%
% OUTPUT;
% X: computed solution to (**)
% solves: number of solves performed to reach the desired tolerance
% res: convergence hystory in terms of residual norm
% err: computed error
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



% check the input
if iscell(A) || iscell(C) || iscell(B)
    error('A, C and B must all be defined as matrices')
elseif ~iscell(N) || ~iscell(M)
    error('N and M must all be defined as cells')
end

% Compute the Schur decomposition once for all
[QA,RA] = schur(full(A));
if norm(size(A)-size(B'))==0 && norm(A-B','fro')<1e-12
        QB=QA;
        RB=RA';
else
    [QB,RB] = schur(full(B));
end

% first solve
Y=lyap(RA,RB,-QA'*C*QB);
    
m = length(N);
nA = size(A,1);
nB = size(B,1);

Ytemp = spalloc(nA,nB,nnz(Y));
NN = cell(m);
MM = NN;
for i=1:m
    % compute once for all QA'*N_i*QA, QB'*M_i*QB
    NN{i}=QA'*N{i}*QA;
    MM{i}=QB'*M{i}*QB;
    
    % compute the product N_i * Y_0 * M_i^T
    Ytemp = Ytemp+NN{i}*Y*MM{i}';
end

solves = 1;    
normC=norm(C,'fro');

res(solves) = norm(Ytemp,'fro')/normC;

while res(solves)>tol && solves<max_solves
    
    % next solve
    YY=lyap(RA,RB,Ytemp);
    
    % update the "reduced" solution
    Y=Y+YY;
    
    % update the number of solves
    solves = solves+1;
    
    % compute the product N_i * Y_j * M_i^T
    Ytemp = spalloc(nA,nB,nnz(Y));
    for i=1:m
        Ytemp = Ytemp+NN{i}*YY*MM{i}';
    end
   
    % compute the residual norm
    res(solves) = norm(Ytemp,'fro')/normC;

end

% compute the actual solution
X=QA*Y*QB';

fprintf('Computed relative residual norm: %10.5e \n',res(end))

% compute the error if the exact solution is available
if nargin>7
    err = norm(Xexact-X,'fro');
    fprintf('Computed error: %10.5e \n',err)  
end
