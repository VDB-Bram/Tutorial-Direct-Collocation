function [C, Ceq] = mycon_DC(z,t1,t2,N,x0,params)

% Optimalisatie Parameters
Tb = z(1);
% d  = z(2);
B  = z(2);
% l  = z(3);
% m  = z(4);
% q  = z(5:4+N);
% qd = z(5+N:4+N*2);
q  = z(3:2+N);
qd = z(3+N:2+N*2);

% Input
deltaT  = (t2-t1)/N;
I = params.I; 
g = params.g;
m = params.m;
l = params.l;
% I       = m*l*l;
% g       = 9.81;
% qdd     = (- m*g*l*cos(q(1:end)))/I + Tb/I + -d*[x0(2) qd(1:end-1)]/I;
qdd     = (- m*g*l*cos(q(1:end)))/I + Tb/I + -B*[x0(2) qd(1:end-1)]/I;

% params.knee_r_range = [-110 0]*pi/180;
% knee_r_range = params.knee_r_range;
% klim = 3;
% theta_ref = (knee_r_range(2) - knee_r_range(1))/2;
% j = x0(1) - (knee_r_range(1) + knee_r_range(2))/2;
% Tlimit = 0*exp(klim*(j-theta_ref)) + 0*exp(klim*(-j-theta_ref));

% OUTPUT 
C    = 0;
Ceq1 = (q-[x0(1) q(1:N-1)])/deltaT - [x0(2) qd(1:end-1)]; 
Ceq2 = qdd - ((qd-[x0(2) qd(1:end-1)])/deltaT);
Ceq  = [Ceq1 Ceq2]';

% Ceq2 = (- mass*g*lc*cos([x0(1) q(1:end-1)]))/I + Tlimit/I + Tb/I - d*[x0(2) qd(1:end-1)]/I - ((qd-[x0(2) qd(1:end-1)])/deltaT);

end 