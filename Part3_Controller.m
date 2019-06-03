%% Part 3 : Adding a controller 
clear all
close all
clc

% Input
data  = load('DataPendulum.mat');
q_exp = data.data(:,2)*pi/180;
t_exp = data.data(:,1);

qdot_exp = diff(q_exp)./diff(t_exp);
qdot_exp(end+1) = qdot_exp(end);

m  = 2.3351;
lc  = 0.2367;
gv = 9.81; 
RG = lc*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc;

% Torques
Tb  = 2;
B   = 0.2;    
kr  = 0.75;
kdr = 0.025;

% Initial state
q0    = q_exp(1);
qdot0 = 0;
q_init = [q0 qdot0]';

% Time info
dt       = 0.0001;                  % time step for backward euler
t_end    = 10;                      % final simulation time
t_vect   = 0:dt:t_end;              % vector with all the time steps
t_span   = t_vect;

params.m  = m;
params.g  = gv;
params.lc = lc;
params.I  = I; 
params.Tb = Tb;
params.B  = B;
params.kr = kr;
params.kdr= kdr;
 
options = odeset('InitialStep',0.01,'MaxStep',0.01);

% Ode
[tM,qM] = ode23(@qdotfunctie_controller, t_span, q_init,options, params,q_exp, qdot_exp);

% ode from part 2 (no controller)
[tM2,qM2] = ode23(@qdotfunctie_discretization, t_span, q_init,options, params);

figure()
plot(tM,qM(:,1))
hold on
plot(tM2,qM2(:,1),'--k')
xlabel('Time [s]'); 
ylabel('Angle [rad]');
legend({'Q: simulated with controller','Q: simulated no controller'})
