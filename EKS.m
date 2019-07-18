classdef EKS < handle
% EKS - Extended Krylov subspace.
%   V=EKS(A,B,m,varargin) progressively builds a basis of an extended
%   Krylov subspace and updates quantities necessary for solving Lyapunov
%   equations, obtaining Ritz vectors, evaluating functions of a matrix,
%   and so on.
%
%   A can be either a matrix or a function handle returning the action of
%   the matrix on a vector. In the latter case, the Ainverse and isSym 
%   flags should be passed to varargin.
%
%   Options for varargin:
%       Ainverse    - allows user to pass in a function handle that applies
%			the inverse of A to a vector.
%       isSym       - Denotes whether this matrix is symmetric.
%       fullOrth    - Fully orthogonalize the basis (in the symmetric case).
%
%   Useful methods:
%       growBasis   - Grow basis in both directions.
%       Vm          - Basis of EKm(A,B).
%       Tm          - Rayleigh quotient matrix Vm'*(A*Vm).
%   For remaining methods, see comments in the code.

    properties
        A;          % input matrix
        Ainverse;   % input matrix
        Ainv;       % for storing invert A
        n;          % size of the space
        r;          % number of columns of B (block size)
        V;          % columns form a basis
        R;          % for Vm'*B
        W;          % storage for new basis vectors to be orthogonalized
        m;          % maximum number of blocks
        j;          % current iterate
        H;          % orthogonalization coefficients
        T;          % projection of A on Krylov subspace
        isSym;      % boolean stating whether A is symmetric or not
        fullOrth;   % do we want to fully orthogonalize?
    end
    
    methods
        function obj=EKS(A,B,m,varargin)
            obj.fullOrth=false;
			for j=1:2:length(varargin)
                switch varargin{j}
                    case 'fullOrth'
                        obj.fullOrth=varargin{j+1};
                    case 'Ainverse'
                        obj.Ainverse=varargin{j+1};
                    case 'isSym'
                        obj.isSym=varargin{j+1};
                end
			end
			if isnumeric(A)
                obj.A = @(X) A*X;
            else
                obj.A = A;
			end
			if isempty(obj.Ainverse)
                obj.Ainv = invertA(A);
                obj.Ainverse = @(X) obj.Ainv.apply(X);
			end
            obj.n = size(A,1);
            obj.m = m;
            obj.r = size(B,2);
            obj.V = zeros(obj.n,2*obj.r*(m+1));
            obj.H = zeros(2*obj.r*(m+1),2*obj.r*m);
            obj.T = zeros(2*obj.r*(m+1),2*obj.r*m);
            obj.W = zeros(obj.n,2*obj.r);
            obj.j = 0;
            [obj.V(:,1:2*obj.r),obj.R] = qr(full([B obj.Ainverse(B)]),0);
            if isempty(obj.isSym)
                obj.isSym = norm(A-A','fro') < 1e-15;
            end
        end
        
        function growBasis(this)
            if this.j==this.m, error('Maximum size exceeded!'), end
            j=this.j+1; r=this.r;
            indsjm1=2*r*(j-2)+1:2*r*(j-1);
            indsj=2*r*(j-1)+1:2*r*j;
            indsjp1=2*r*j+1:2*r*(j+1);
			
			% Lanczos
			if this.isSym && ~this.fullOrth
                inds=@(j,k) max(1,2*r*(j-2)+k);
			else
                inds=@(j,k) 1;
			end
			
			% Grow block as in Simoncini 2007
            this.W=[this.A(this.V(:,indsj(1:r))) this.Ainverse(this.V(:,indsj(r+1:end)))];
			for k=1:2*r
                for reorth=1:2
                    for i=inds(j,k):2*r*j+k-1
                        h=this.W(:,k)'*this.V(:,i);
                        this.H(i,2*r*(j-1)+k)=this.H(i,2*r*(j-1)+k)+h;
                        this.W(:,k)=this.W(:,k)-h*this.V(:,i);
                    end
                end
                this.H(2*r*j+k,2*r*(j-1)+k)=norm(this.W(:,k));
                this.V(:,2*r*j+k)=this.W(:,k)/this.H(2*r*j+k,2*r*(j-1)+k);
			end
            
			% Update Rayleigh quotient matrix T = V'*A*V as in Simoncini 2007
            if j==1
                this.T(1:4*r,1:r)=this.H(1:4*r,1:r);
                this.T(1:4*r,r+1:2*r)=( [this.R(1:r,1:r);zeros(3*r,r)] - ...
                    this.H(1:4*r,1:r)*this.R(1:r,r+1:2*r) ) / this.R(r+1:2*r,r+1:2*r);
            else
                if this.isSym && ~this.fullOrth
                    Ipiece=kron([0;1;0;0;0;0],eye(r));
                else
                    sn=2*(j+1); sk=2*(j-1);
                    Ipiece=kron( [zeros(sk-1,1);1;zeros(sn-sk,1)] ,eye(r));
                end
                this.T(inds(j,1):indsjp1(end),indsj(1:r))=this.H(inds(j,1):indsjp1(end),indsj(1:r));    
                this.T(inds(j,1):indsjp1(end),indsj(r+1:end))=( Ipiece - ...
                    this.T(inds(j,1):indsjp1(end),1:indsj(r))*this.H(1:indsj(r),indsjm1(r+1:end)) ) ...
                    / this.H(indsj(r+1:end),indsjm1(r+1:end));
            end
            
            this.j = this.j+1;
        end
        
        function Vm=Vm(this)
            Vm=this.V(:,1:this.j*2*this.r);
        end
        
        function Vmp1=Vmp1(this)
            Vmp1=this.V(:,1:(this.j+1)*2*this.r);
        end
        
        function vmp1=vmp1(this)
            vmp1=this.V(:,2*this.r*this.j+1:(this.j+1)*2*this.r);
        end
        
        function Hm=Hm(this)
            Hm=this.H(1:this.j*2*this.r,1:this.j*2*this.r);
        end
        
        function Hmbar=Hmbar(this)
            Hmbar=this.H(1:(this.j+1)*2*this.r,1:this.j*2*this.r);
        end
        
        function Tm=Tm(this)
            Tm=this.T(1:this.j*2*this.r,1:this.j*2*this.r);
        end
        
        function Tmbar=Tmbar(this)
            Tmbar=this.T(1:(this.j+1)*2*this.r,1:this.j*2*this.r);
        end
        
        function tmp1m=tmp1m(this)
            tmp1m=this.T(2*this.r*this.j+1:(this.j+1)*2*this.r,2*this.r*(this.j-1)+1:this.j*2*this.r);
        end
        
        function x=blockSize(this)
            x=2*this.r;
        end
    end
end