function [qdot] = qdotfunction_controller (t_span,q_init, params, q_exp, qdot_exp)

% Input 
m  = params.m;
g  = params.g;
lc = params.lc;
I  = params.I; 
Tb = params.Tb;
B  = params.B;
kr = params.kr;
kdr= params.kdr;

% Dynamics 
qdot    = zeros(2,1);
qdot(1) = q_init(2);
qdot(2) = (-m*g*lc*cos(q_init(1)))/I + Tb/I + -B*q_init(2)/I + (kr*(q_init(1)-q_exp(1))+kdr*(q_init(2)-qdot_exp(1)))/I;

