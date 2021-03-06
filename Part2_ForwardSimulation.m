%% Part 2 : Discretization - Backward Euler %%
% This script contains information about a forward simulations (see
% chapter2 on forward simulations) of the pendulum test. Three different
% approaches are introduced.
%
%   2A: time-marching with a backward euler integration scheme
%   2B: time-marching using matlab Ordinary Differential Equation solvers (ODE)
%   2C: Discretise equations and solve numerically using fsolve in matlab

% Notes about symbols:

% q      joint agle
% qdot   angular velocity
% qddot  angular acceleration
% x      state: consists of joint angle and angular velocity


%% Path information

% add functions to the matlab path
addpath(fullfile(pwd,'functions'));

%% Load pendulum properties (see Part1_Input.m for more information)

% Input
data  = load('DataPendulum.mat');
q_exp = data.data(:,2)*pi/180;
t_exp = data.data(:,1);

m  = 2.3351;            % mass of the pendulum
lc = 0.2367;            % distance knee - COM tibia
gv = 9.81;              % gravitational constant
RG = lc*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc; % Inertia of lower limb +foot

% Torques
Tb = 2;
B  = 0.2;              % B will be used to calculate Td(damping) Td=-B*qdot

% Initial state of the pendulum
x0 = [q_exp(1); 0];  % [joint angle initial frame  -  zero velocity ]


%% Part 2A : Time marching (for loop)

% Time steps
dt       = 0.005;                    % time step for backward euler
t_end    = 10;                      % final simulation time
t_vect   = 0:dt:t_end;              % vector with all the time steps
N        = length(t_vect);          % number of steps

% Pre allocate states
x = nan(2,length(t_vect));          % pre-allocate matrix with states

% set the initial state
x(:,1) = x0;

for i=1:N-1    
    % get angles and velocities at selected time point
    q    =  x(1,i);
    qdot = x(2,i);    
    % equations of motion of the pendulum (see Chapter 1: Dynamics of the pendulum test)
    T         = Tb + B*qdot;                    % Torques
    qddot     = (-m*gv*lc*cos(q) - T)./I  ;     % angular velocity    
    % Backward euler:
    %  q(t+dt) - q(t) 
    %  -------------    = qdot(t)
    %        dt    
    qi        = qdot*dt+q;
    qdoti     = qddot*dt +qdot;    
    % store the predicted state
    x(1,i+1) = qi;
    x(2,i+1) = qdoti;
    
end

% plot the predicted angle and angular velocity
figure();
lw = 2;     % linewidth
subplot(1,2,1);
plot(t_vect,x(1,:),'LineWidth',lw);xlabel('Time [s]'); ylabel('Angle [rad]');
hold on 
title('Forward simulation of knee angle')
subplot(1,2,2);
plot(t_vect,x(2,:),'LineWidth',lw); xlabel('Time [s]'); ylabel('Angular velocity [rad/s]');
hold on 
legend({'time marching'})
title('Forward simulation of knee angular velocity')

%% part 2B: Ode 

% Time info
t_span   = [0 t_end];    % time interval for ode

params.m  = m;
params.g  = gv;
params.lc = lc;
params.I  = I; 
params.Tb = Tb;
params.B  = B;
 
% some options for the ode soler
options =  odeset('RelTol',1e-8,'AbsTol',1e-8);

% using ODE function in matlab to solve Ordinary Differential Equation.
% Note that we created the function "StateDerivative" to compute the time
% derivative of x. Open this function to gain more insight in how this
% works.
[t,x] = ode23(@StateDerivative, t_span, x0 ,options, params);

% add this to the figure
subplot(1,2,1);
plot(t,x(:,1),'g','LineWidth',lw)

subplot(1,2,2);
plot(t,x(:,2),'g','LineWidth',lw)



%% Part 2C: Implicit (fsolve)

%  store the model parameters in the structure "params"
params.m  = m;
params.lc = lc;
params.g  = gv;
params.I  = I;
params.Tb = Tb;
params.B  = B;

% time vector
dt       = 0.005;                  % time step for backward euler
t_end    = 10;                      % final simulation time
t_vect   = 0:dt:t_end;  

% pre-allocate states
x = zeros(2,length(t_vect));

% solve set of equations using fsolve matlab
% Note: Open the function "mycon_fsolve" to gain insight in how we construct the
% set of equations, wich will be solved suing the fsolve function.
F_handle = @(x)mycon_fsolve(x,dt,x0,params);
x        = fsolve(F_handle,x);
q        = x(1,:);
qdot     = x(2,:);

subplot(1,2,1);
plot(t_vect,q,'--r','LineWidth',lw)

subplot(1,2,2);
plot(t_vect,qdot,'--r','LineWidth',lw)
legend({'time marching','ode23','fsolve'})


%% Additional notes

% You can gain some insight on the error you make on the integration with
% the backward euler scheme depending on the time step. Change the time
% step dt to 0.01 for example and evaluate how this influences your
% results. You can use the solution with the ODE solver as golden standard
% since the integration error here is minimal due to the variable step
% (Note that we set the error tolerance low for the ODE solve = 1e-8).

