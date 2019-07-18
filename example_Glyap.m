function example
n=70;
h=1/(n+1);
e=ones(n,1);
e1=sbv(n,1);
en=sbv(n,n);
E1=sparse(e1*e1');
En=sparse(en*en');
K=spdiags([-e 2*e -e],-1:1,n,n);
C=spdiags([-e 0*e e],-1:1,n,n);
I=speye(n);

A=( - kron(I,K) - kron(K,I) + .5*( kron(E1,I) + kron(En,I) ) )/h^2 + ...
	kron(I,C)/(2*h);
N=.5*[ kron(E1,I) kron(En,I) ]/h;
B=-[ kron(e1,e) kron(en,e) ]/h;

	tol=struct('inexact',1e-2,'outer',1e-8);
	maxit=struct('inner',50,'outer',20);
	method.res=1;
	method.rhs=1;
	[Z,data]=efficientGenLyap(A,N,B,tol,maxit,method);


function e=sbv(n,k)
%SBV   Standard basis vector.
%   E = SBV(N,K) returns the Kth standard unit basis vector of dimension N,
%   i.e., a column vector of size N of zeros with a 1 in the Kth entry.

e=zeros(n,1);
e(k)=1;

