function control = mi20_control()
% function control = mi20_control()
%
% The control structure is initilized. The components directly relate to
% the control structure in the HSL_MI20 documentation.
%
% After initialization, the components may by altered. For example, 
% if 10 v-cycle iterations are required, use the following command:
%    control.v_iterations = 10;

control = hsl_mi20('control');