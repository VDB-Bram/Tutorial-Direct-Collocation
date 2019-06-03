function [Ceq] = mycon_fsolve(x,dt,q_init,params)
q = x(:,1);
qdot = x(:,2);

Tb = params.Tb;
B  = params.B;
m  = params.m;
lc = params.lc;
g  = params.g;
I  = params.I;

% Dynamics
qdd     = (- m*g*lc*cos(q(1:end)))/I + Tb/I + -B*[q_init(2); qdot(1:end-1)]/I;

% OUTPUT 
Ceq1 = (q-[q_init(1); q(1:end-1)])/dt - [q_init(2); qdot(1:end-1)]; 
Ceq2 = qdd - ((qdot-[q_init(2); qdot(1:end-1)])/dt);
Ceq  = [Ceq1; Ceq2]';

end

