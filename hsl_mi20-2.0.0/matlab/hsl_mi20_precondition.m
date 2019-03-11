function [x,inform] = hsl_mi20_precondition(z,i)
% function [x,inform] = hsl_mi20_precondition(z)
%
% Having previously called mi20_setup, the AMG preconditioner is 
% applied to the vector z and x is returned.
%
% The function optionally returns a struture inform that contains 
% information from the execution of mi20_precondition. 
%  If inform.flag==0, the code successfully executed with no
%    warnings or errors.
%  If inform.flag>0, the code successfully executed but a warning
%    was created. 
%  If inform.flag<0, then and error occured.
%  The value of inform.flag corresponds to the value of the
%  info%flag in the HSL_MI20 documentation.

if nargin == 1
   i = 1; 
end

global mi20_handle;
global MI20STATUS;
 
if (MI20STATUS(i) ~= 1)
  error('mi20_setup must be called first')
end 

[x,inform] =  hsl_mi20( 'precondition',z,mi20_handle{i});

