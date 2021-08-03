function plotmesh( x, z )
%PLOTMESH Plots a 2D grid of mesh lines and mesh line density.

%   Ben Postlethwaite 2012
%   benpostlethwaite.ca

x = x(:);
z = z(:);
x = x';
z = z';

nz = length(z);
nx = length(x);

line([x; x], [ones(nx,1)' * min(z) ; ones(nx,1)' *max(z)],'Color',[.8 .8 .8])
line([ones(nz,1)' * min(x) ; ones(nz,1)' *max(x)], [z; z], 'Color',[.8 .8 .8])

xlim( [min(x) - 0.01*min(x), max(x) + 0.01*min(x)] )
ylim( [min(z) - 0.01*min(z), max(z) + 0.01*min(z) ] )
title('2D Mesh')
set(gca,...
'XTickLabel','', 'XMinorTick', 'off')
set(gca,'XTick',[])
axis square
end

