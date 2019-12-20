function [xdot] = StateDerivative(t_span,x,params)

% Input 
m  = params.m;
g  = params.g;
lc = params.lc;
I  = params.I; 
Tb = params.Tb;
B  = params.B;

% get angles and velocities
q = x(1);
qdot =x(2);

% Dynamics 
xdot    = zeros(2,1);
xdot(1) = qdot;
xdot(2) = (-m*g*lc*cos(q))/I - Tb/I -B*qdot/I ;

