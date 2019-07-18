function [Z,data]=efficientGenLyap(A,N,B,tol,maxit,method)
% EFFICIENTGENLYAP   EKSM based solver for generalized Lyapunov equations.
%	[Z,data]=EFFICIENTGENLYAP(A,N,B,tol,maxit,method) solves a generalized
%	Lyapunov equation (or GLE)
%		A*X + X*A' + N_1*X*N_1' + ... + N_m*X*N_m' + B*B' = 0
%	numerically for a low-rank factor Z such that X \approx Z*Z'. The input
%	N must be of the form N = [N_1 , ... , N_m].
%
%	Input tol should be a struct. The field outer is the stopping criteria
%	for the GLE. We let rho_k denote an estimate of the outer residual norm
%	which may be exact (see below). If tol.inexact is set, we will an use
%	an inner or truncation tolerance of the form tol.inexact*rho_k,
%	depending on whether or not these tolerances are specified. If
%	tol.inner or tol.trunc is set, we used a fixed inner or truncation
%	tolerance.
%
%	Input maxit should be a struct with fields inner and outer.
%
%	The struct method determines how the inner Lyapunov solve is done and
%   how the residual is checked. rho_k is is the true residual norm if
%	method.res=0, or a bound (described in SSS14) if method.res=1. A
%	"thick" right-hand side is used if method.rhs=0, and a separating of
%	the right-hand described (described in SS14) is used if method.rhs=1.
%
%	Output is a low-rank factor Z of an approximation to X and a struct
%	data with fields relres, solves, memory, rankZ, and time.

tic
% Initialization steps
data=struct('relres',[],'solves',0,'memory',[],'rankZ',[]);
normBBT=normGDGT(B);			% Calculate || B*B' ||_F inexpensively
ZtruncTol=1e-14;				% Tolerance for progressive truncation of Z
rB=size(B,2);					% rank of B
if method.res, rNZold=0; end	% Rank of stored previous Z; nothing stored initially
rZ=0;							% for memory computation at initial inner iterate of initial outer iterate
m=size(N,2)/size(A,1);			% number of N_j terms
Ainv=invertA(A);				% factorize A
Ainverse=@(X) Ainv.apply(X);	% application of A^{-1} to vector; pass to EKS

for k=1:maxit.outer     % Outer iteration
    % Calculate rhs for next EKSM solve, set tolerances
    if k==1
        % Inner tolerance for first solve
		if isfield(tol,'inner')
			innerTol=tol.inner;
		else
			innerTol=tol.inexact;
		end
		if isfield(tol,'trunc')
			truncTol=tol.trunc;
		else
			truncTol=tol.inexact;
		end
        rhs=B;    % First rhs for extended Krylov subspace
    else
        % Solves represents total number of solves; record from previous iterates
        data.solves(k)=data.solves(k-1);
        
        % Inner tolerance for next solve
		if isfield(tol,'inner')
			innerTol=tol.inner;
		else
			innerTol=tol.inexact*data.relres(k-1);
		end
        
        % Truncation tolerance for rhs
		if isfield(tol,'trunc')
			truncTol=tol.trunc;
		else
			truncTol=tol.inexact*data.relres(k-1);
		end

        % Truncated low-rank factor B_k of B_k*B_k' := N_1*X_k*N_1^T + ... + N_m*X_k*N_m^T + B*B'
        [rhs,rhsNorm]=truncationZZT(full([NN(N,Z) B]),truncTol);
    end
    normBkBkT=normGDGT(rhs);	% Norm of rhs
    rRhs=size(rhs,2);			% Rank of rhs
    
    % Inner iteration; Lyapunov solve via EKSM
	fprintf('\tINNER SOLVE -- Rank of rhs: %d\n',rRhs)
	if method.res	% Initial memory demands, depending on method
		stage1mem=rNZold;
	else
		stage1mem=0;
	end

	if method.rhs==1
		% Solve each Lyapunov equation invidually
		fprintf('\tRhs ')
		for i=1:rRhs
			fprintf('%d ',i)
			% Initialize extended Krylov subspace
			V=EKS(A,rhs(:,i),maxit.inner,'Ainverse',Ainverse);
			beta=V.R(1,1);                      % For Vm'*B*B'*Vm = beta^2*e1*e1'
			data.solves(k)=data.solves(k)+1;    % One solve done in above
			for j=1:maxit.inner
				V.growBasis                         % Expand basis of Krylov subspace
				data.solves(k)=data.solves(k)+1;    % One solve done in above

				% Solve projected equation
				Xtilde=lyap(V.Tm,diag([beta^2;zeros(2*j-1,1)]));

				% Calculate residual norm cheaply
				absres=sqrt(2)*norm(V.tmp1m*Xtilde(2*(j-1)+1:2*j,:),'fro');
				relres=absres/(rRhs*normBkBkT);
				if relres < innerTol/rRhs   % Check for termination (oversolving)
					break
				end
			end

			% Assemble Z
			if i==1
				% Xtilde should be symmetric; symmetrize for accuracy
				[Q,D]=eig((Xtilde+Xtilde')/2);

				% Discard negligible parts
				ind=find(abs(diag(D))/max(abs(diag(D)))>ZtruncTol);
				ZQ=V.Vm*Q(:,ind);
				ZR=sqrt(D(ind,ind));

				% Update memory usage
				if method.res
					stage1mem=rNZold+size(V.Vm,2);
				else
					stage1mem=size(V.Vm,2);
				end
			else
				% Xtilde should be symmetric; symmetrize for accuracy
				[Q,D]=eig((Xtilde+Xtilde')/2);

				% Discard negligible parts
				ind=find(abs(diag(D))/max(abs(diag(D)))>ZtruncTol);
				VmQ=V.Vm*Q(:,ind);
				Y=D(ind,ind);

				% Update memory usage
				if method.res
					stage1mem=max(stage1mem,rNZold+rZ+size(V.Vm,2));
				else
					stage1mem=max(stage1mem,rZ+size(V.Vm,2));
				end

				% Progressively truncate
				ZQTV=ZQ'*VmQ;
				[Q2,R2]=qr(VmQ-ZQ*ZQTV,0);
				core=[ZR*ZR'+ZQTV*Y*ZQTV'	ZQTV*Y*R2'
					  R2*Y*ZQTV'			R2*Y*R2'];
				% core should be symmetric; symmetrize for accuracy
				[Q,D]=eig((core+core')/2);

				% Discard negligible parts
				if i<rRhs
					ind=find(abs(diag(D))/max(abs(diag(D)))>ZtruncTol);
				else
					[~,ind1]=sort(abs(abs(diag(D))),'descend');
					Q=Q(:,ind1);
					D=D(ind1,ind1);
					maxD= max(abs(diag(D)));
					ind2=find(abs(diag(D))/maxD> truncTol,1,'last');
					ind=1:ind2;
				end

				ZQ=[ZQ Q2]*Q(:,ind);
				ZR=sqrt(D(ind,ind));
				clear Q2
				rZ=size(ZQ,2);
			end
			clear V
		end
		Z=ZQ*ZR;
		rZ=size(ZQ,2);
		clear ZQ ZR
	else
		% Solve with a thick rhs
		% Initialize extended Krylov subspace
		V=EKS(A,rhs,maxit.inner,'Ainverse',Ainverse);
		beta=V.R(1:rRhs,1:rRhs);				% For Vm'*B*B'*Vm = beta*E1*E1'*beta'
		data.solves(k)=data.solves(k)+rRhs;		% rRhs solves done in above
		fprintf('Inner it:')
		for j=1:maxit.inner
			fprintf(' %d',j)
			V.growBasis							% Expand basis of Krylov subspace
			data.solves(k)=data.solves(k)+rRhs;	% rRhs solves done in above

			% Solve projected equation
			Xtilde=lyap(V.Tm,kron(sbv(2*j,1)*sbv(2*j,1)',beta*beta'));

			% Calculate residual norm cheaply
			absres=sqrt(2)*norm(V.tmp1m*Xtilde(2*rRhs*(j-1)+1:2*rRhs*j,:),'fro');
			relres=absres/normBkBkT;
			if relres < innerTol	% Check for termination (oversolving)
				break
			end
		end
		% Xtilde should be symmetric; symmetrize for accuracy
		[Q,D]=eig((Xtilde+Xtilde')/2);

		% Discard negligible parts
		ind=find(abs(diag(D))/max(abs(diag(D)))>truncTol);
		Z=V.Vm*Q(:,ind)*sqrt(D(ind,ind));
		rZ=size(Z,2);	% Rank of latest Z
	end
	fprintf('\n')
    
    % Memory demands (part 2) and residual norms
    if method.res
		if k==1
            NZ=NN(N,Z);
			NZnorm=normGDGT(NZ);
            data.relres(k)=(innerTol*normBBT+NZnorm)/normBBT;
            stage2mem=rZ+size(NZ,2);
        else
            NZmZold=[NN(N,Z) NZold];
			NZmZoldnorm=normGDGT(NZmZold,blkdiag(eye(rZ*m),-eye(rNZold)));
            data.relres(k)=(truncTol*rhsNorm+innerTol*normBkBkT+ ...
                NZmZoldnorm)/normBBT;
            stage2mem=rZ+size(NZmZold,2);
		end
    else
        stage2mem=rZ+3*rZ+rB;
        data.relres(k)=normGDGT([Z A*Z NN(N,Z) B], ...
            blkdiag( kron( blkdiag([0 1;1 0],eye(m)) , eye(rZ) ) , eye(rB) ) ) / normBBT;
    end
    data.memory(k)=max(stage1mem,stage2mem);    % Larger of two memory demands
    data.rankZ(k)=rZ;                           % Rank of low-rank factor
    
    % Display progress, stop if converged
    fprintf('\nOUTER ITERATE %d:\trelres:%d\trank:%d\tsolves:%d\n', ...
        k,data.relres(k),data.rankZ(k),data.solves(k))
    fprintf('  Memory:%d\tStage 1:%d\tStage 2:%d\tRank rhs:%d\n\n',data.memory(k), ...
        stage1mem,stage2mem,rRhs)
    if data.relres(k) < tol.outer, break, end
    
    % Keep low-rank factor for residual checking (if needed)
	if method.res
        NZold=NN(N,Z);
        NZold=truncationZZT(NZold,ZtruncTol);
        rNZold=size(NZold,2);
	end
end
data.time=toc;
fprintf('TOTAL CPU TIME:%.2f\n',data.time)

function x=normGDGT(G,D)
% NORMGDGT   Norm of a low-rank symmetric matrix.
%	Given low-rank factors G and D of a symmetric matrix X = G*D*G',
%	return || X ||_F in an efficient manner.

if nargin<2, D=eye(size(G,2)); end
[~,R]=qr(full(G),0);
x=norm(R*D*R','fro');

function Y=NN(N,Z)
% NN   Apply operator \mathcal{N}, described in SSS14.
%	Given a low-rank factor Z of a matrix X = Z*Z', return low-rank factor
%	Y of Y*Y' = N_1*X*N_1' + ... + N_m*X*N_m'.

n=size(N,1);	% size of each N_j
m=size(N,2)/n;	% number of N_j terms
r=size(Z,2);	% number of columns of low-rank factor Z
Y=zeros(n,m*r);
for j=1:m
	Y(:,1+(j-1)*r:j*r)=N(:,1+(j-1)*n:j*n)*Z;
end

function [Z,normZ0] = truncationZZT(Z0,truncTol)
% TRUNCATIONZZT   Truncation of a low-rank symmetric positive semidefinite
%   matrix.
%   Z=truncationZZT(Z0,TRUNCTOL) returns a matrix Z such that
%       || Z*Z' - Z0*Z0' ||_F <= TRUNCTOL*|| Z0*Z0' ||_F,
%   where the rank of Z0 is the smallest integer such that the above
%   inequality holds. Also returns || Z0*Z0' ||_F.

[Q,R]=qr(Z0,0);
RRT=R*R';
normZ0=norm(RRT,'fro');
[V,D]=eig((RRT+RRT')/2);
S=abs(diag(D));
[~,ind]=sort(S,'ascend');
norms=[0;sqrt(cumsum(S(ind).^2))/normZ0];
k=find(norms<truncTol,1,'last');
Z=Q*V(:,ind(k:end))*sqrt(D(ind(k:end),ind(k:end)));


function e=sbv(n,k)
%SBV   Standard basis vector.
%   E = SBV(N,K) returns the Kth standard unit basis vector of dimension N,
%   i.e., a column vector of size N of zeros with a 1 in the Kth entry.

e=zeros(n,1);
e(k)=1;

