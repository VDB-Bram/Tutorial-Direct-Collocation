%% OptimalisatieParametersPendulumJente
clear all
clc
close all

%% Notes
% Optimizes the damping and baseline torque of an inverted pendulum model
% to track human data.

%% Load data
% Experimental Data
data  = load('DataPendulum.mat');
t_exp = data.data(:,1);
q_exp = data.data(:,2);

% Input parameters
m  = 2.318;
lc = 0.2368;
g  = 9.81;
RG = lc*0.416;   % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc;

%% Load Casadi and linear solver
%addpath('C:\software\linear_solver');
% Casadi = 'C:\software\Casadi'; addpath(genpath(Casadi));
import casadi.*;

%% Time information for mesh- and collocation points

dt     = 0.01; % time between mesh points
t_mesh = t_exp(1):dt:t_exp(end)+dt;
N      = length(t_mesh)-1; % number of mesh points
t_coll = linspace(t_mesh(1),t_mesh(end),N*3);


qspline              = spline(t_exp,q_exp);
% q_exp at meshpoints
[qmeshExp, qdmeshExp, qddmeshExp] = SplineEval_ppuval(qspline,t_mesh,1);
% q_exp at collocationpoints
[qcollExp, qdcollExp] = SplineEval_ppuval(qspline,t_coll,1);

%% Ipopt options

options.ipopt.tol = 1*10^(-6);
options.ipopt.linear_solver = 'mumps';
options.ipopt.hessian_approximation = 'limited-memory';

%% Casadi

% setup polynomials
d = 3;                  % degree of polynomials
tau_root = [0 collocation_points(d, 'radau')];
C = zeros(d+1,d+1);     % Coefficients of the collocation equation
D = zeros(d+1, 1);      % Coefficients of the continuity equation
B = zeros(d+1, 1);      % Coefficients of the quadrature function

% Construct polynomial basis
for j=1:d+1
    % Construct Lagrange polynomials to get the polynomial basis at the collocation point
    coeff = 1;
    for r=1:d+1
        if r ~= j
            coeff = conv(coeff, [1, -tau_root(r)]);    % multiplication of the polynomials
            coeff = coeff / (tau_root(j)-tau_root(r));
        end
    end
    % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D(j) = polyval(coeff, 1.0);
    
    % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the collocation equation
    pder = polyder(coeff);
    for r=1:d+1
        C(j,r) = polyval(pder, tau_root(r));
    end
    
    % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
    pint = polyint(coeff);
    B(j) = polyval(pint, 1.0);
end

%% Opti
opti = casadi.Opti();

% Bounds
Bounds.q        = [-5 5];
Bounds.qdot     = [-30 30];
Bounds.qddot    = [-100 100];

% Mesh points
q    = opti.variable(1,N+1);    % joint angles
qd   = opti.variable(1,N+1);    % angular velocity
qdd  = opti.variable(1,N);      % acceleration

opti.subject_to(Bounds.q(1) < q < Bounds.q(2));
opti.subject_to(Bounds.qdot(1) < qd < Bounds.qdot(2));
opti.subject_to(Bounds.qddot(1) < qdd < Bounds.qddot(2));

% Collocation points
qcoll       = opti.variable(1,N*d);
qdcoll      = opti.variable(1,N*d);

opti.subject_to(Bounds.q(1) < qcoll < Bounds.q(2));
opti.subject_to(Bounds.qdot(1) < qdcoll < Bounds.qdot(2));

% Initial guess
% Meshpoints
opti.set_initial(q,qmeshExp);
opti.set_initial(qd,qdmeshExp);
opti.set_initial(qdd,qddmeshExp(1:end-1));

% Collocationpoints
opti.set_initial(qcoll,qcollExp);
opti.set_initial(qdcoll,qdcollExp);

% Design variables
% Bounds
Bounds.uDamp    = [0.01 1];
Bounds.uTb      = [-10 10];

uDamp           = opti.variable(1,1);       % Damping coefficient
opti.subject_to(Bounds.uDamp(1) < uDamp < Bounds.uDamp(2));

uTb           = opti.variable(1,1);         % BaselinTorque
opti.subject_to(Bounds.uTb(1) < uTb < Bounds.uTb(2));

% Initial guess
opti.set_initial(uDamp,0);
opti.set_initial(uTb,0);

% Bounds on initial state
opti.subject_to(q(1) - q_exp(1) == 0);      % set initial angles
opti.subject_to(qd(1) == 0);             % set initial velocities to zero

% Formulate the NLP
% loop over meshpoints
for k=1:N
    % get the state and controls
    qs       = q(k);
    qds      = qd(k);
    qs_coll  = [qs qcoll((k-1)*d+1:k*d)];
    qds_coll = [qds qdcoll((k-1)*d+1:k*d)];
    qdds = qdd(k);
    
    % loop over the collocation points
    for j=1:d
        % Expression of the state derivatives at the collocation points
        qp   = qs_coll*C(:,j+1);
        qdp  = qds_coll*C(:,j+1);
        
        opti.subject_to(dt*qds_coll(:,j+1) - qp == 0);
        opti.subject_to(dt*qdds  - qdp == 0);
    end
    
    % Skeleton dynamics
    opti.subject_to(qdds == (-m*g*lc*cos(qs)- uDamp*qds + uTb)./I);
    
    % state continuity at mesh transition
    opti.subject_to(q(k+1) - qs_coll*D == 0);
    opti.subject_to(qd(k+1)- qds_coll*D == 0);
end

% objective function
qerror = q-qmeshExp;
J      = sumsqr(qerror);
opti.minimize(J);

% Solve
opti.solver('ipopt',options);
sol = opti.solve();

% Solutions
qsol        = sol.value(q)';
qdsol       = sol.value(qd)';
uDampsol    = sol.value(uDamp);
uTbsol      = sol.value(uTb);

%% Plot results

figure()
plot(t_mesh,qsol)
hold on
plot(t_mesh,qmeshExp,'--k')
hold on
legend({'q sol','q exp'})
