%% Direct Collocation - fmincon
clear all
close all
clc

%% Load data

% Experimental Data
data  = load('DataPendulum.mat');
t_exp = data.data(:,1);
q_exp = data.data(:,2)*pi/180;

f1    = 1;     % first frame of simulation
f2    = 101;   % final frame of simulation
t_exp = t_exp(f1:f2,:);
q_exp = q_exp(f1:f2,:);

t1       = t_exp(1); 
t2       = t_exp(101);  
tspan    = [t1 t2]';
N        = length(0:0.01:(t2-t1));
t        = t_exp(1);

% Input
m  = 2.3351;
l  = 0.2367;
g  = 9.81; 
RG = l*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*l*l;

% Initial State
x0 = [q_exp(1) 0]'; % initial state [q qdot TSRS Tact]

% Qdot for initial guess   
qdot_exp        = diff(q_exp)/0.01; 
qdot_exp(end+1) = qdot_exp(end);
    
%% Params
params.N     = N; 
params.q_exp = q_exp;
params.m     = m;
params.l     = l;
params.g     = g;
params.I     = I;
params.t     = t; 

%% DC
% [z, fval ] = fmincon(F_hanlde,ig,A,B,Aeq,Beq,Lb,Ub,nonlcon);
    % z = parameters geoptimaliseerd
    % fval = uitkomst van minimialisatie
    % F_hanlde = functie die geminimaliseerd moet worden (difference
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
F_hanlde2 = @(z)mycon_DC(z,t1,t2,N,x0,params);
options   = optimoptions('fmincon','MaxFunctionEvaluations',1000000);

% Initial Guess
load('C:\Users\u0125183\Documents\MATLAB\Tutorial Numerieke Simulatie/qMinterp.mat')
ig(1)         = 0;                  % ig(Tb)
ig(2)         = 0;                  % ig(B)
ig(3:2+N)     = qMinterp(1:101,1);  % ig(q)
ig(3+N:2+2*N) = qMinterp(1:101,2);  % ig(qdot)

% Equality constraints 
A = [];     Aeq= [];
B = [];     Beq= [];

% Bounds 
Lb(1)       = 0;            Ub(1)       = 10;
Lb(2)       = 0;            Ub(2)       = 10;
Lb(3:2+2*N) = -300;         Ub(3:2+2*N) = 300;

% Simulation 
[z, fval ] = fmincon(F_hanlde,ig,A,B,Aeq,Beq,Lb,Ub,F_hanlde2,options);

%% Plot

figure()
plot(t_exp,z(3:2+N))
hold on
plot(t_exp,q_exp, '--k')
legend({'Q; DC with Fmincon','Q : Experimental'})
