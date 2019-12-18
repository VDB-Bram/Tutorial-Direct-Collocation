function [Ceq] = mycon_fsolve(x,dt,q_init,params)

% angles
qt1 = x(:,1);                           % t+dt
q   = [q_init(1); qt1(1:end-1)];        % t

% velocities
qdott1 = x(:,2);                        % t+dt          
qdot = [q_init(2); qdott1(1:end-1)];    % t

% parameters pendulum model
Tb = params.Tb;
B  = params.B;
m  = params.m;
lc = params.lc;
g  = params.g;
I  = params.I;

% Dynamics
qdd     = (- m*g*lc*cos(q))/I - Tb/I + -B*qdot/I;

% OUTPUT (equality constraints
Ceq1 = (qt1-q)/dt - qdot; 
Ceq2 = qdd - ((qdott1-qdot)/dt);
Ceq  = [Ceq1; Ceq2]';

end

