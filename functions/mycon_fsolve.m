
function [Ceq] = mycon_fsolve(x,dt,x0,params)

% Input arguments:
% x: states (q and qdot) at discretized time
% dt: time step
% x0: state at time 0
% params: structure with model properties

% Output:
% Ceq: vector with constraint equations
N = length(x(1,:));

% parameters pendulum model
Tb = params.Tb;
B  = params.B;
m  = params.m;
lc = params.lc;
g  = params.g;
I  = params.I;

%% Vector implementation:
% Note: this is harder to understand, but computationally more efficient

% angles
qt1 = x(1,2:N);         % angles at time: t+dt
q   = x(1,1:N-1);       % angle at time: t

% velocities
qdott1 = x(2,2:N);         % angles at time: t+dt
qdot   = x(2,1:N-1);       % angle at time: t

% Dynamics of the pendulum
qdd     = (- m*g*lc*cos(q))/I - Tb/I + -B*qdot/I;

% constraints (integration scheme)
Ceq1 = qdot - (qt1-q)/dt;           % backward euler velocity 
Ceq2 = qdd - ((qdott1-qdot)/dt);    % backward euler acceleration
Ceq  = [Ceq1; Ceq2];               % vector with constraints

% constraints on initial state
Cx0 = x(:,1) - x0;

% append constrains
Ceq = [Ceq Cx0];


%% For loop implementatation (easier to understand)
% Note: this is easier to understand, but computationally less efficient

% % pre allocate vector with constraints
% Ceq1 = zeros(1,N-1);
% Ceq2 = zeros(1,N-1);
% 
% % for loop
% for i=1:N-1
%     % get angles and velocities at time t
%     q       = x(1,i);
%     qdot    = x(2,i);
%     
%      % get angles and velocities at time t+1
%     qt1       = x(1,i+1);
%     qdott1    = x(2,i+1);
%     
%     % Dynamics of the pendulum
%     qdd     = (- m*g*lc*cos(q))/I - Tb/I + -B*qdot/I;
% 
%     % OUTPUT (equality constraints
%     Ceq1(i) = qdot - (qt1-q)/dt;           % backward euler velocity 
%     Ceq2(i) = qdd - ((qdott1-qdot)/dt);    % backward euler acceleration    
% end
% Ceq  = [Ceq1; Ceq2]';               % vector with constraints
% % constraints on initial state
% Cx0 = x(:,1) - x0;
% 
% % append constrains
% Ceq = [Ceq Cx0];

end

