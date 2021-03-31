function func = meshfunc(curve, weight)
%MESHFUNC Returns function of set parmaters for building graded meshes.
%   Meshfunc takes a set of parameters and an optional interactive mode
%   flag and returns a function for building graded meshes. The returned
%   function takes a regularly spaced vector x and returns a two vectors.
%
%   gmesh = meshfunc(curve, wieght, interactive)
%   [x, dx] = gmesh( linspace(0, 1, n) )
%
%   x samples values on the graded mesh and dx is equal to diff(x)
%
%   The Inputs to meshfunc are curve - which adjusts the steepness in the
%   change of dx spacing. Weight, which adjusts the ratio between the
%   larger dx spacing at the boundaries and the smallest dx spacing in the
%   centre dx(centre) / dx(boundary).

%   Ben Postlethwaite 2012
%   benpostlethwaite.ca

c = curve;
w = weight;

% Gmesh builds a fat parabola with minimum dictated by wieght
% The w/(1-w) is necessary to that the scaling is preserved.
gmesh = @(x) ((x(:) - (x(end) - x(1))/2 - x(1)).^c + (w / (1-w) ) * ((x(end) - x(1))/2).^c);
% getdx turns gmesh into an x and dx
getdx = @(x) (x(end) - x(1)) * gmesh(x) / sum(gmesh(x));
%getdx = @(x) gmesh(x);
getx = @(x) cumsum([x(1); getdx(x)]);
% Project cell walls to centres
centre = @(x) x(1:end - 1) + 0.5*diff(x);
% Get dims performs a few function calls to make things easier.
func = @(x) deal( getx(x), getdx(x) ); % Removed centre(x) from both getx and getdx calls

end



