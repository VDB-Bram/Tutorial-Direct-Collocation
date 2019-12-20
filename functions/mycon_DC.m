function [C, Ceq] = mycon_DC(z,dt,params)

N = params.N;

% Optimalisation Parameters
Tb  = z(1);
B   = z(2);
q   = z(3:2+N);
qdot = z(3+N:2+N*2);

% Input
I       = params.I; 
g       = params.g;
m       = params.m;
lc      = params.lc;

% Dynamics
qddot     = (-m*g*lc*cos(q))/I - Tb/I -B*qdot/I ;

% equality constraints based on backward euler scheme
Ceq1 = (q(2:N) - q(1:N-1))./dt - qdot(1:N-1);
Ceq2 = (qdot(2:N) - qdot(1:N-1))./dt - qddot(1:N-1);

% output equality constraints
C    = [];
Ceq  = [Ceq1 Ceq2]';

%------------------------------------------------
% NOTE: It is easier to understand the integration scheme using a for loop
% notation
% Ceq1    = zeros(N-1,1);
% Ceq2    = zeros(N-1,1);
% for i=1:N-1
%     Ceq1(i) = (q(i+1)-q)./dt - qdot(i);
%     Ceq2(i) = (qdot(i+1)-qdot)./dt - qddot(i);
% end


end 