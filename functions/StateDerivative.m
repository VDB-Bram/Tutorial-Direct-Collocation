function [xdot] = StateDerivative(t,x,params)
%StateDerivative compuates the time derivative of the state x 
% input: 
%   (1) t: time
%   (2) x: state (consists of [q qdot])
%   (3) params: structure with parameters of the pendulum

% output:
%   (1) xdot: time derivative of x ([qdot qddot])


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

