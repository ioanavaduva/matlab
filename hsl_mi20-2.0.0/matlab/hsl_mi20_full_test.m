function hsl_mi20_full_test()
% Test code for hsl_mi20 Matlab interface

disp( 'Checking correct errors' )

A = gallery('wathen',30,30);
b = ones(length(A),1);

% Error test 1
% No inputs
try
    x = hsl_mi20();
catch
    errstr = lasterror;
    str=strtrim(errstr.message);
    str1=strtrim('Insufficient input arguments');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 1')
    end
end

% Error test 2
% Input string doesn't match entry from allowable list of strings
try
    x = hsl_mi20('cnt');
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Unrecognised action: ''cnt''');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 2')
    end
end

% Error test 3
% Too many output arguments
try
    [x,y] = hsl_mi20( 'control' );
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('hsl_mi20 provides at most 1 output argument for this call');
    if (size(strfind(str,str1),2)==0)
       errstr
       errstr.message
        error('Failure at error test 3')
    end
end

% Error test 4
% Not enough input arguments for 'setup'
try
    [x,y] = hsl_mi20( 'setup' );
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Insufficient input arguments');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 4')
    end
end

% Error test 5
% Not enough input arguments for 'setup'
try
    [x,y] = hsl_mi20( 'setup',A );
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Insufficient input arguments');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 5')
    end
end

% Error test 6
% Too many input arguments for 'setup'
try
    [x,y] = hsl_mi20( 'setup',A,1,1 );
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('argument ''control'' should be a structure');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 6')
    end
end

% Error test 7
% Not enough input arguments for 'precondition'
try
    [x,y] = hsl_mi20( 'precondition', 1 );
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Insufficient input arguments');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 7')
    end
end


% Error test 8
% Invalid handle
try
    [x,y] = hsl_mi20( 'precondition',1,-1 );
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Invalid handle');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 8')
    end
end

% Error test 9
% wrong type of input argument
try
    [x,y] = hsl_mi20( A );
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('First argument must be string');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 9')
    end
end

control = hsl_mi20('control');
control.coarse_solver=3;
[i,handle] = hsl_mi20( 'setup', A,control);


% Error test 10
% too many input arguments for 'precondition'
try
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Wrong number of input arguments to hsl_mi20');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 10')
    end
end

% Error test 11
% too many output arguments
try
    [x,i,y] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Too many output arguments');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 11')
    end
end


% Error test 12
% second input argument must be a vector
try
    [x,i] =  hsl_mi20( 'precondition',ones(length(A),2),handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Error in argument z. Expected 1-dimensional matrix.');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 12')
    end
end


% Error test 13
%second input argument must be a vector
try
    [x,i] =  hsl_mi20( 'precondition',sparse(ones(length(A),2)),handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('z must only have one column');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 13')
    end
end

% Error test 14
% cd must be a
try
    [x,inform] = hsl_mi20( 'precondition',ones(length(A),1),handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('There must be a vector data');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 14')
    end
end


% Error test 15
% A must have positive diagonals
try
    A(1,1)=0;
    [i,handle] = hsl_mi20( 'setup',A,control);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('One or more diagonal entry missing');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 15')
    end
end


% Error test 16
% A must have positive diagonals
try
    A(1,1)=-1;
    [i,handle] = hsl_mi20( 'setup',A,control);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('One or more diagonal entry is <= 0');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 16')
    end
end

A = gallery('wathen',30,30);

% Error test 17
% control.st_parameter too large
try
    control.st_parameter = 1.1;
    [i,handle] = hsl_mi20( 'setup',A,control);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.st_parameter out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 17')
    end
end

control = hsl_mi20('control');
% Error test 18
% control.err_tol negative
try
    control.err_tol = -1.1;
    [i,handle] = hsl_mi20( 'setup',A,control);
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.err_tol out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 18')
    end
end

control = hsl_mi20('control');
% Error test 19
% control.err_tol negative
try
    control.max_points = -1.1;
    [i,handle] = hsl_mi20( 'setup',A,control);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.max_points out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 19')
    end
end

control = hsl_mi20('control');
% Error test 20
% control.err_tol negative
try
    control.st_method = -1.1;
    [i,handle] = hsl_mi20( 'setup',A,control);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.st_method out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 20')
    end
end


control = hsl_mi20('control');
% Error test 21
% control.aggressive negative
try
    control.aggressive = -1.1;
    [i,handle] = hsl_mi20( 'setup',A,control);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.aggressive out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 21')
    end
end

control = hsl_mi20('control');
% Error test 22
% control.c_fail out of range
try
    control.c_fail = 4;
    [i,handle] = hsl_mi20( 'setup',A,control);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.c_fail out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 22')
    end
end

control = hsl_mi20('control');
% Error test 23
% control.smoother zero
try
    control.smoother=0;
    [i,handle] = hsl_mi20( 'setup',A,control);
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.smoother out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 23')
    end
end

control = hsl_mi20('control');
% Error test 24
% control.pre-smoothing negative
try
    control.pre_smoothing=-1;
    [i,handle] = hsl_mi20( 'setup',A,control);
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.pre_smoothing out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 24')
    end
end

control = hsl_mi20('control');
% Error test 25
% control.post-smoothing negative
try
    control.post_smoothing=-1;
    [i,handle] = hsl_mi20( 'setup', A,control);
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.post_smoothing out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 25')
    end
end

control = hsl_mi20('control');
% Error test 26
% control.pre_smoothing+control.post_smoothing = 0
try
    control.pre_smoothing=0;
    control.post_smoothing=0;
    [i,handle] = hsl_mi20( 'setup',A,control);
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.pre_smoothing+control.post_smoothing = 0');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 26')
    end
end

control = hsl_mi20('control');
% Error test 27
% control.coarse_solver=0
try
    control.coarse_solver=0;
    [i,handle] = hsl_mi20( 'setup', A,control);
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.coarse_solver out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 27')
    end
end

control = hsl_mi20('control');
% Error test 28
% control.coarse_solver_its negative
try
    control.coarse_solver=0;control.coarse_solver_its = -1;
    [i,handle] = hsl_mi20( 'setup', A,control);
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.coarse_solver_its out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 28')
    end
end

control = hsl_mi20('control');
% Error test 29
% control.damping out of range
try
    control.smoother = 1;control.damping = -1;
    [i,handle] = hsl_mi20( 'setup',A,control);
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.damping out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 29')
    end
end

control = hsl_mi20('control');
% Error test 30
% control.max_levels out of range
try
    control.max_levels = -1;
    [i,handle] = hsl_mi20( 'setup', A,control);
    [x,i] =  hsl_mi20( 'precondition',b,handle);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('control.max_levels out of range');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 30')
    end
end

display('Finished checking errors')
clear
%hsl_mi20_startup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
('Testing one application of AMG applied to the Poisson problem:\n') ;

A = gallery('poisson',20,20);
b = A*ones(length(A),1);
control = hsl_mi20_control;
inform = hsl_mi20_setup(A,control);
x =  hsl_mi20_precondition(b);

fprintf ('error in 2-norm: %f\n \n',norm(x-A\b)) ;

hsl_mi20_finalize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hsl_mi20_startup;

fprintf ('Testing AMG in conjunction with PCG applied to the Poisson problem:\n') ;
control = hsl_mi20_control;

inform = hsl_mi20_setup(A,control);
x = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
hsl_mi20_finalize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Correct execution test 3

%hsl_mi20_startup;

control = hsl_mi20_control;
control.aggressive = 2;
inform = hsl_mi20_setup(A,control);
[x,flag] = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
if (flag~=0)
    error('Error in correct execution test 3')
end
hsl_mi20_finalize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Correct execution test 4

%hsl_mi20_startup;

control = hsl_mi20_control;
control.c_fail = 2;
inform = hsl_mi20_setup(A,control);
[x,flag] = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
if (flag~=0)
    error('Error in correct execution test 4')
end
hsl_mi20_finalize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Correct execution test 5

A = gallery('poisson',30,30);
b = A*ones(length(A),1);
%hsl_mi20_startup;

control = hsl_mi20_control;
control.max_levels = 20;
inform = hsl_mi20_setup(A,control);
[x,flag] = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
if (flag~=0)
    error('Error in correct execution test 5')
end
hsl_mi20_finalize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Correct execution test 6

%hsl_mi20_startup;

control = hsl_mi20_control;
control.max_points = 5;
inform = hsl_mi20_setup(A,control);
[x,flag] = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
if (flag~=0)
    error('Error in correct execution test 6')
end
hsl_mi20_finalize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Correct execution test 7

%hsl_mi20_startup;

control = hsl_mi20_control;
control.reduction = 0.5;
inform = hsl_mi20_setup(A,control);
[x,flag] = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
if (flag~=0)
    error('Error in correct execution test 7')
end
hsl_mi20_finalize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Correct execution test 8

%hsl_mi20_startup;

control = hsl_mi20_control;
control.st_method = 1;
inform = hsl_mi20_setup(A,control);
[x,flag] = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
if (flag~=0)
    error('Error in correct execution test 8')
end
hsl_mi20_finalize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Correct execution test 9

%hsl_mi20_startup;

control = hsl_mi20_control;
control.coarse_solver = 2;
inform = hsl_mi20_setup(A,control);
[x,flag] = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
if (flag~=0)
    error('Error in correct execution test 9')
end
hsl_mi20_finalize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Correct execution test 10

%hsl_mi20_startup;

control = hsl_mi20_control;
control.pre_smoothing = 4;
inform = hsl_mi20_setup(A,control);
[x,flag] = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
if (flag~=0)
    error('Error in correct execution test 10')
end
hsl_mi20_finalize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Correct execution test 11

%hsl_mi20_startup;

control = hsl_mi20_control;
control.v_iterations = 4;
inform = hsl_mi20_setup(A,control);
[x,flag] = pcg(A,b,1e-8,100,'hsl_mi20_precondition');
if (flag~=0)
    error('Error in correct execution test 11')
end
hsl_mi20_finalize;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Test for more than one instance of AMG 13
clear all
n_vec = [20,30];

for i = 1:length(n_vec)
    n = n_vec(i);
    A{i} = gallery('poisson',n,n);
    b{i} = A{i}*ones(length(A{i}),1);
    control{i} = hsl_mi20_control;
    inform = hsl_mi20_setup(A{i},control{i},i);
    pre{i} = @(x) hsl_mi20_precondition(x,i);
end

for i = 1:length(n_vec)
    [x{i},flag] = pcg(A{i},b{i},1e-8,100,pre{i});
    if (flag~=0)
        error('Error in more than one instance test 13')
    end
end

hsl_mi20_finalize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Test with non-sequential preconditioners 14
n_vec = [20,30];

for i = 1:length(n_vec)
    in = 2*(i-1) + 1;
    n = n_vec(i);
    A{i} = gallery('poisson',n,n);
    b{i} = A{i}*ones(length(A{i}),1);
    control{i} = hsl_mi20_control;
    inform = hsl_mi20_setup(A{i},control{i},in);
    pre{i} = @(x) hsl_mi20_precondition(x,in);
end

for i = 1:length(n_vec)
    [x{i},flag] = pcg(A{i},b{i},1e-8,100,pre{i});
    if (flag~=0)
        error('Error in more than one instance test 15')
    end
    hsl_mi20_finalize(2*(i-1) + 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test with individual finalizes 15

n_vec = [20,30];

for i = 1:length(n_vec)
    n = n_vec(i);
    A{i} = gallery('poisson',n,n);
    b{i} = A{i}*ones(length(A{i}),1);
    control{i} = hsl_mi20_control;
    inform = hsl_mi20_setup(A{i},control{i},i);
    pre{i} = @(x) hsl_mi20_precondition(x,i);
end

for i = 1:length(n_vec)
    [x{i},flag] = pcg(A{i},b{i},1e-8,100,pre{i});
    if (flag~=0)
        error('Error in more than one instance test 15')
    end
    hsl_mi20_finalize(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Test with finalizing a non-existing handle 16
%hsl_mi20_finalize(2);
n_vec = [20,30];

for i = 1:length(n_vec)
    in = 2*(i-1) + 1;
    n = n_vec(i);
    A{i} = gallery('poisson',n,n);
    b{i} = A{i}*ones(length(A{i}),1);
    control{i} = hsl_mi20_control;
    inform = hsl_mi20_setup(A{i},control{i},in);
end

% test finalizing a non-initialized array
try
    hsl_mi20_finalize(2);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('The index sent to the function does not correspond to an instance of hsl_mi20');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 16')
    end
end

hsl_mi20_finalize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Test with finalizing a too high handle 16

try
    hsl_mi20_finalize(2);
catch
    errstr=lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('The index sent to the function does not correspond to an instance of hsl_mi20');
    if (size(strfind(str,str1),2)==0)
        error('Failure at error test 16')
    end
end


display('All tests succeeded.')