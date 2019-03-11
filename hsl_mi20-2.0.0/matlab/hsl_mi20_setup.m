function inform = hsl_mi20_setup(A,control,i)
% function inform = hsl_mi20_setup(A,control[,i])
%
% Given a matrix A, mi20_setup will setup the AMG preconditioner
%
% control is a structural argument that is initialized with the call
% control = mi20_control()
%
% the optional parameter i is an integer that identifies an instance of
% the amg preconditioner, allowing more than one to be used in the
% code at a time
%
% The function returns a struture inform that contains information from 
% the execution of mi20_setup. 
%  If inform.flag==0, the code successfully executed with no warnings or errors
%  If inform.flag>0, the code successfully executed but a warning was created. 
%  If inform.flag<0, then an error occured.
%   The value of i.flag corresponds to the value of the info%flag in the 
%     HSL_MI20 documentation.
%
% Subsequent calls to mi20_setup may only be made after mi20_finalize has
% been called
if nargin == 2
    i = 1;
end

global mi20_handle; 
global MI20STATUS;


[inform,mi20_handle{i}] = hsl_mi20( 'setup', A, control);

MI20STATUS(i) = 1;
