%% Optimization using shooting and Fmincon 

clear all
close all
clc

% Input
data  = load('DataPendulum.mat');
t_exp = data.data(:,1);
q_exp = data.data(:,2)*pi/180;
t_span = [t_exp(1) t_exp(end)]';

m  = 2.3351;
lc = 0.2367;
gv = 9.81; 
RG = lc*0.416;               % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc;

% Params
global params 
params.m     = m;
params.lc    = lc;
params.g     = gv;
params.I     = I;

% Initial state
q_init = [q_exp(1) 0]';      %[q0 qdot0]

% Initial guess
ig = [0 0.01];               % Tb B 
lb = [0 0.01] ;              % Tb B 
ub = [10 10] ;               % Tb B  
    
%% Optimization
    
% Optimization of Tb,B to minimize the difference between predicted
% and experimental values 

% Function that should be minimized = myobj_shooting 
F_handle  = @(z)myobj_shooting(z,q_exp,t_span,q_init,t_exp',params);
[z, fval] = fmincon(F_handle,ig,[],[],[],[],lb,ub);

%% Validation
    
params.Tb = z(1);
params.B  = z(2);
options = odeset('InitialStep',0.01,'MaxStep',0.01);

[tM,qM] = ode23(@qdotfunctie_shooting, t_span, q_init,options, params);
   
figure();
plot(tM, qM(:,1)); hold on;
plot(t_exp,q_exp,'--k');
legend({'Q: Shooting','Q: Experimental'})
title('Shooting with optimized Tb and B ')
xlabel('Time [s]'); 
ylabel('Angle [rad]');
