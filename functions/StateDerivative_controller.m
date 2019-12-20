function [xdot] = StateDerivative_controller(t,x, params, qdes, qdot_des)

% get angles and velocities
q   = x(1);
qd  = x(2);

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
xdot    = zeros(2,1);
xdot(1) = x(2);

% torques:
Tdamp = B*qd;
Tfb   = kr*(q-qdes)+kdr*(qd-qdot_des); 
Ttot  = -Tb - Tdamp - Tfb;

% get the joint acceleration
qdd  = (-m*g*lc*cos(q))/I + Ttot/I;

% append state
xdot(2) = qdd;

