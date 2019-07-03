%RKFUNB    RKFUNB constructor.
%
% Usages:
% 1) obj = rkfunb(K, H, coeffs, k) constructs an RKFUNB from the
%    recursion matrices (H,K), and the coefficient matrices 
%    provided in coeffs.
%    The optional parameter k specifies sub- or superdiagonal
%    rational functions.



classdef rkfunb
properties
K
H
coeffs
k
end % properties
    
methods(Access = public, Static = true)
% Constructor. (in methods)
function obj = rkfunb(varargin)
  
  % Convert from RKFUN.
  if isa(varargin{1},'rkfun')
    obj.K = varargin{1}.K;
    obj.H = varargin{1}.H;
    obj.coeffs = num2cell(varargin{1}.coeffs);
    obj.k = varargin{1}.k;
    return
  end
  
  % Convert from RKFUNM.
  if isa(varargin{1},'rkfunm')
    [m,n] = size(varargin{1}.coeffs(1));
    obj.K = kron(varargin{1}.K,eye(m));
    obj.H = kron(varargin{1}.H,eye(m));
    obj.coeffs = varargin{1}.coeffs;
    obj.k = varargin{1}.k;
    return
  end
  
    
  % Gallery.
  if ischar(varargin{1})
    obj = rkfunb.gallery(varargin{:});
    return
  end
  
  % Construct rkfunb from pencil and coeffs.
  obj.K = varargin{1};
  obj.H = varargin{2};
  obj.coeffs = varargin{3};
  if nargin >= 4
    obj.k = varargin{4};
  else
    obj.k = 0;
  end
  
end % function rkfunb
        
function obj = gallery(varargin)
% GALLERY    Collection of rational matrix-valued functions.
%
% obj = rkfunm.gallery(funname, param1, param2, ...) takes 
% funname, a case-insensitive string that is the name of 
% a rational function family, and the family's input 
% parameters. 
%
% At the moment the gallery is empty. 

  obj = gallery(varargin{:});

end

end % methods    
end % classdef rkfun
