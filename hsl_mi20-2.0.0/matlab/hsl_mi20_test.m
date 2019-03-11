%hsl_mi20_startup

fprintf ('Testing one application of AMG applied to the Poisson problem:\n') ;

A = gallery('poisson',20,20);
b = A*ones(length(A),1);
control = hsl_mi20_control;
inform = hsl_mi20_setup(A,control);
x = hsl_mi20_precondition(b);

fprintf ('error in 2-norm: %f\n \n',norm(x-A\b)) ;
x = hsl_mi20_precondition(x);
x = hsl_mi20_precondition(x);
hsl_mi20_finalize;

fprintf ('\n\nTesting AMG in conjunction with PCG applied to the Poisson problem:\n') ;
control = hsl_mi20_control;
inform = hsl_mi20_setup(A,control);
x = pcg(A,b,1e-8,100,'hsl_mi20_precondition');

hsl_mi20_finalize;

clear A b control inform x 

fprintf ('\n\nTesting more than one instance of AMG at once\n\n')

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
    x{i} = pcg(A{i},b{i},1e-8,100,pre{i});
end

hsl_mi20_finalize;