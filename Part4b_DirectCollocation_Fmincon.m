%% Direct Collocation - fmincon

%% add functions to the matlab path
addpath(fullfile(pwd,'functions'));

%% Load pendulum properties (see Part1_Input.m for more information)

% Input
data  = load('DataPendulum.mat');
q_exp = data.data(:,2)*pi/180;
t_exp = data.data(:,1);

% t_span = [t_exp(1) t_exp(end)]';
% note use smaller t_span, because this implementation is very slow 
t_span = [t_exp(1) t_exp(1)+2]';    % simulate 2s

% pendulum properties
m  = 2.3351;            % mass of the pendulum
lc = 0.2367;            % distance knee - COM tibia
gv = 9.81;              % gravitational constant
RG = lc*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc; % Inertia of lower limb +foot

% Params
params.m     = m;
params.lc    = lc;
params.g     = gv;
params.I     = I;

% Initial state of the pendulum
x0 = [q_exp(1); 0];  % [joint angle initial frame  -  zero velocity ]

%% DC
% [z, fval ] = fmincon(F_hanlde,ig,A,B,Aeq,Beq,Lb,Ub,nonlcon);
    % z = optimization parameters
    % fval = output of the optimization
    % F_handle = objective function (i.e. difference between measured and
    % simulation joint angles)
    % ig = initial guess of all paramaters
    % A = Linear equality constraints (vector)
    % B = Linear inequality constraints (vector)
    % Aeq = Linear equality constraints (matrix)
    % Beq = Linear equality constraints (matrix)
    % Lb = Lower bound of optimized parameters
    % Ub = Upper bound of optimized parameters
    % nonlcon = Non linear constraints 
    
% Notes
%  (1) The vector z contains all the optimization variables. The first two
% elements are the static parameters Tb and B. After this, the joint angles
% on the discretized times and the angular velocities are stored.
% (2) We selected the 

% Initial guess damping and baseline torque
ig(1)         = 0;     % Initial guess Tb
ig(2)         = 0;     % Initial guess B

% discretisation for collocation
dt      = 0.01;                    % time step
tvect   = t_span(1):dt:t_span(2);  % time vector (discretisation)
N       = length(tvect);

% interpolate experimenatal data at discr. time using spline
qspline = spline(t_exp,q_exp);  % spline fit
[q_exp,qdot_exp] = SplineEval_ppuval(qspline,tvect,1); % get angles and velocities

% set initial guess
ig(3:2+N)     = q_exp;  % Initial guess joint angle
ig(3+N:2+2*N) = qdot_exp;  % Initial gues angular velocity

% linear constraints (empty because we used nonlinear constraint
A = [];     Aeq= [];
B = [];     Beq= [];

% Bounds on the optimization variables
Lb = zeros(size(ig));       Ub = zeros(size(ig));
Lb(1)       = -10;          Ub(1)       = 10;   % baseline torqu between 0 and 10
Lb(2)       = -10;          Ub(2)       = 10;   % damping between 0.01 and 10
Lb(3:end)   = -300;         Ub(3:end)   = 300;    %  angles and angular velocities between -300 and 300 

% constraint on initial state (fixed bound)
Lb(3)   = x0(1);   Lb(3+N)   = x0(2);       % Upper bound on angles and angular velocity respectively
Ub(3)   = x0(1);   Ub(3+N)   = x0(2);       % Lower bound on angles and angular velocity respectively

% add additional information to params
params.N = N;
params.q_exp = q_exp;

% objective function and constraints equations
f_obj       = @(z)myobj_DC(z,params);       % objective function 
f_constr    = @(z)mycon_DC(z,dt,params); % constraints function

% Solve optimization problem using fmincon
options     = optimoptions('fmincon','MaxFunctionEvaluations',1000000);
[z, fval ]  = fmincon(f_obj,ig,A,B,Aeq,Beq,Lb,Ub,f_constr,options);

%% Plot

figure()
plot(tvect,z(3:2+N))
hold on
plot(tvect,q_exp, '--k')
legend({'Q: DC with Fmincon','Q: Experimental'})
xlabel('Time [s]'); 
ylabel('Angle [rad]');
title('DC - fmincon: optimized Tb and B ')
