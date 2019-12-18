function [qdot] = StateDerivative(t_span,q_init, params)

% Input 
m  = params.m;
g  = params.g;
lc  = params.lc;
I  = params.I; 
Tb = params.Tb;
B  = params.B;

% Dynamics 
qdot    = zeros(2,1);
qdot(1) = q_init(2);
qdot(2) = (-m*g*lc*cos(q_init(1)))/I - Tb/I + -B*q_init(2)/I ;

