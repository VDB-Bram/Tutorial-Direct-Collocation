%% Optimalisation using shooting and Fmincon 

clear all
close all
clc

%% Load data
% Experimental Data
data  = load('DataPendulum.mat');
t_exp = data.data(:,1);
q_exp = data.data(:,2)*pi/180;
t_span = [t_exp(1) t_exp(end)]';

% Input
m  = 2.3351;
l  = 0.2367;
g  = 9.81; 
RG = l*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*l*l;

% Input to params
global params 
params.m     = m;
params.lc    = l;
params.g     = g;
params.I     = I;

% Initial state
q_init = [q_exp(1) 0]';

% Initial guess
ig = [0 0];       %Tb B 
lb = [0 0.01] ;   %Tb B 
ub = [10 10] ;    %Tb B  
    
%% Optimization
    
% Optimization of Tb,B to minimize the difference between predicted
% and experimental values 

% Function that should be minimized = myobj 
F_hanlde  = @(z)myobj_shooting(z,q_exp,t_span,q_init,t_exp',params);
[z, fval] = fmincon(F_hanlde,ig,[],[],[],[],lb,ub);

%% Validation
    
params.Tb = z(1);
params.B  = z(2);
options = [];

[tM,qM] = ode23(@qdotfunctie_shooting, t_span, q_init,options, params);
   
figure();
plot(tM, qM(:,1)); hold on;
plot(t_exp,q_exp,'--k');
legend({'Shooting','experimental'})
title('Shooting with optimized Tb and B ')

%% Save output

qMinterp = (interp1(tM,qM,t_exp));
save('C:\Users\u0125183\Documents\MATLAB\Tutorial Numerieke Simulatie/qMinterp.mat','qMinterp')
