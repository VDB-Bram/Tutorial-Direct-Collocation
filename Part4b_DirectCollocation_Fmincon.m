%% Direct Collocation - fmincon
clear all
close all
clc

% Input
data  = load('DataPendulum.mat');
t_exp = data.data(:,1);
q_exp = data.data(:,2)*pi/180;

% f1    = 1;     % first frame of simulation
% f2    = 501;   % final frame of simulation
% t_exp = t_exp(f1:f2,:);
% q_exp = q_exp(f1:f2,:);

f1 = 1
f2 = length(t_exp)

t1       = t_exp(f1); 
t2       = t_exp(f2);  
tspan    = [t1 t2]';
N        = length(0:0.01:(t2-t1));
t        = t_exp(1);

m  = 2.3351;
lc = 0.2367;
gv = 9.81; 
RG = lc*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc;

% Initial State
q_init = [q_exp(1) 0]'; % initial state [q0 qdot0]

% Qdot for initial guess   
qdot_exp        = diff(q_exp)/0.01; 
qdot_exp(end+1) = qdot_exp(end);
    
%% Params
params.N     = N; 
params.q_exp = q_exp;
params.m     = m;
params.lc    = lc;
params.g     = gv;
params.I     = I;
params.t     = t; 

%% DC
% [z, fval ] = fmincon(F_hanlde,ig,A,B,Aeq,Beq,Lb,Ub,nonlcon);
    % z = parameters geoptimaliseerd
    % fval = uitkomst van minimialisatie
    % F_handle = functie die geminimaliseerd moet worden (difference
    % between q and theta(i) ^2
    % ig = initial guess of all paramaters
    % A = Linear equality constraints (vector)
    % B = Linear inequality constraints (vector)
    % Aeq = Linear equality constraints (matrix)
    % Beq = Linear equality constraints (matrix)
    % Lb = Lower bound of optimized parameters
    % Ub = Upper bound of optimized parameters
    % nonlcon = Non linear constraints 

F_hanlde  = @(z)myobj_DC(z,params,t);
F_hanlde2 = @(z)mycon_DC(z,t1,t2,N,q_init,params);
options   = optimoptions('fmincon','MaxFunctionEvaluations',1000000);

% Initial Guess
    % Simulation without optimizing Tb and B as ig for q and qdot
    params.Tb = 2;
    params.B  = 0.2;
    options_ode = odeset('InitialStep',0.01,'MaxStep',0.01); 
    [tM,qM] = ode23(@qdotfunctie_discretization, tspan, q_init,options_ode, params);

ig(1)         = 0;        % ig(Tb)
ig(2)         = 0.01;     % ig(B)
ig(3:2+N)     = qM(:,1);  % ig(q)
ig(3+N:2+2*N) = qM(:,2);  % ig(qdot)

% Equality constraints 
A = [];     Aeq= [];
B = [];     Beq= [];

% Bounds 
Lb(1)       = 0;            Ub(1)       = 10;
Lb(2)       = 0.01;         Ub(2)       = 10;
Lb(3:2+2*N) = -300;         Ub(3:2+2*N) = 300;

% Simulation 
[z, fval ] = fmincon(F_hanlde,ig,A,B,Aeq,Beq,Lb,Ub,F_hanlde2,options);

%% Plot

figure()
plot(t_exp,z(3:2+N))
hold on
plot(t_exp,q_exp, '--k')
legend({'Q: DC with Fmincon','Q: Experimental'})
xlabel('Time [s]'); 
ylabel('Angle [rad]');
title('DC - fmincon: optimized Tb and B ')
