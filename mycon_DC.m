function [C, Ceq] = mycon_DC(z,t1,t2,N,x0,params)

% Optimalisation Parameters
Tb = z(1);
B  = z(2);
q  = z(3:2+N);
qd = z(3+N:2+N*2);

% Input
deltaT  = (t2-t1)/N;
I       = params.I; 
g       = params.g;
m       = params.m;
lc      = params.lc;

% Dynamics
qdd     = (- m*g*lc*cos(q(1:end)))/I + Tb/I + -B*[x0(2) qd(1:end-1)]/I;

% OUTPUT 
C    = 0;
Ceq1 = (q-[x0(1) q(1:N-1)])/deltaT - [x0(2) qd(1:end-1)]; 
Ceq2 = qdd - ((qd-[x0(2) qd(1:end-1)])/deltaT);
Ceq  = [Ceq1 Ceq2]';

end 