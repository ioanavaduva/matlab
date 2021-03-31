%% Graded Mesh (a dumb octree)
% curve and weight are prechosen after running
% meshgui, then we can just run meshfunc with
% best curve and weight parameters

%   Ben Postlethwaite 2012
%   benpostlethwaite.ca

N = 25;
curve = 4;
weight = 0.1; 
gmesh = meshfunc(curve, weight);

% OR
% To choose curve and weight via the GUI do
% 
% gmesh = meshgui();

% Use gmesh func for creating graded mesh out of 
% linearly spaced inputs
[z, dz] = gmesh(linspace(5, 6, N - 10));
[x, dx] = gmesh(linspace(2, 3, N + 5));
[y, dy] = gmesh(linspace(0, 1, N + 5));

% View mesh in 2D
figure()
subplot(1,2,1)
plotmesh(z, x)

subplot(1,2,2)
plotmesh(x, y)

