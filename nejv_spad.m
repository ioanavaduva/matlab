function [X, res, i] = nejv_spad( A, B, tol, numit)
% [X, res, i] = nejv_spad( A, B, tol, numit)
% input: A - the system matrix (symmetric positive definite, it is not checked)
%        B - the right hand side 
%        tol - tolerance for norm of the residual
%        numit - maximal number of iterations
% output: X the solution of AX=B computed with given tolerance or after numit iterations
%         res - norm of the residual
%         i - number of iterations
  X = zeros(size(B));   % the first approx. (zero vector of size of B)                
  R = B;                % residual
  res = norm(R);        % norm of the residual
  i = 0;
  % iterations
  while ((i < numit) && (res > tol))
    i = i+1;
    Q = A * R';
    alfa = (R * R) / (R' * Q);
    X = X + alfa * R;
    R = R - alfa * Q;
    res = norm(R);
  end  % while
end  % function