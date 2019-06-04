%% OptimalisatieParametersPendulumJente
% script to track experimental data by optimizing baselinetone, damping and
% reflex gains

clear all
clc
%% Input data
% Subject specific data
mtotaal      = 38.28;
lbeen        = sqrt(0.0310*0.0310+0.089*0.089+0.379*0.379);

time_of_drop = 13.88;
tdrop        = time_of_drop;
qinit        = 0.1130;

% Experimental data
q_exp_data   = load('DataPendulum.mat');
q_exp        = q_exp_data.data(:,2)*pi/180;
q_exp_time   = q_exp_data.data(:,1);
m            = 0.061*mtotaal;
l            = 0.606*lbeen;
g            = 9.81;           % gravity
RG           = (l*0.416);
I            = m*RG*RG + m*l*l;

% Time information
t1          = tdrop;
t2          = q_exp_data.data(end,1);
dt          = 0.01;            % time step backward euler
tspan       = [t1 t2]';         % total simulation time
tvect       = 0:dt:(t2-t1);     % vector with the time points
N           = length(tvect);    % total collocation points

q_exp_interp = q_exp(find(q_exp_time==t1):find(q_exp_time==t2));

%% Casadi en linear solver
addpath('C:\software\linear_solver');   % for ma57
Casadi = 'C:\software\Casadi'; addpath(genpath(Casadi));
import casadi.*;

options.ipopt.tol = 1*10^(-6);
options.ipopt.linear_solver = 'ma57';

% copy collocation scheme params to model structure
h     = dt;                % time step
times = tvect;             % time vector

ki_StartObj = 1;
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

knee_r_range = [-110 0]*pi/180;
klim         = 3;
b            = 0.1;

% Bounds
Bounds.uTlim    = [-100 100];
Bounds.q        = [-10 10];
Bounds.qdot     = [-100 100];
Bounds.uDamp    = [0.01 1];
Bounds.uTb      = [-10 10];
Bounds.uM       = [m m];
Bounds.uL       = [l l];
Bounds.uR       = [0.01 1];
Bounds.udotR    = [0.001 1];
Bounds.uRh      = [0.01 1];
Bounds.udotRh   = [0.001 1];

% model variables
q    = opti.variable(1,N+1);
qd   = opti.variable(1,N+1);
qdd  = opti.variable(1,N);

% States collocation points
qcoll       = opti.variable(1,N*d);
qdcoll      = opti.variable(1,N*d);
opti.subject_to(Bounds.q(1) < qcoll < Bounds.q(2));
opti.subject_to(Bounds.qdot(1) < qdcoll < Bounds.qdot(2));

% controls - mesh points
uDamp           = opti.variable(1,1);
opti.subject_to(Bounds.uDamp(1) < uDamp < Bounds.uDamp(2));

uTb           = opti.variable(1,1);
opti.subject_to(Bounds.uTb(1) < uTb < Bounds.uTb(2));

uM            = opti.variable(1,1);
opti.subject_to(Bounds.uM(1) < uM < Bounds.uM(2));

uL            = opti.variable(1,1);
opti.subject_to(Bounds.uL(1) < uL < Bounds.uL(2));

uR            = opti.variable(1,1);
opti.subject_to(Bounds.uR(1) < uR < Bounds.uR(2));

udotR         = opti.variable(1,1);
opti.subject_to(Bounds.udotR(1) < udotR < Bounds.udotR(2));

uRh            = opti.variable(1,1);
opti.subject_to(Bounds.uRh(1) < uRh < Bounds.uRh(2));

udotRh         = opti.variable(1,1);
opti.subject_to(Bounds.udotRh(1) < udotRh < Bounds.udotRh(2));

% Set the initial guess
opti.set_initial(q,zeros(1,N+1));
opti.set_initial(qd,zeros(1,N+1));
opti.set_initial(uDamp,0.1);
opti.set_initial(uTb,zeros(1,1));
opti.set_initial(uM,m);
opti.set_initial(uL,l);
opti.set_initial(uR,0.75);
opti.set_initial(udotR,0.025);
opti.set_initial(uRh,0.75);
opti.set_initial(udotRh,0.025);
opti.set_initial(qcoll,zeros(1,N*d));
opti.set_initial(qdcoll,zeros(1,N*d));

% bounds on initial state
opti.subject_to(q(1)  - qinit == 0);      % set initial angles
opti.subject_to(qd(1) - 0 == 0);          % set initial velocities to zero

J = 0;  % initliase cost function

% Formulate the NLP
for k=1:N
    % get the state and controls
    qs  = q(k);     qds = qd(k);    %uDamps  = uDamp(k); uTbs = uTb(k);
    qs_coll  = [qs qcoll((k-1)*d+1:k*d)];
    qds_coll = [qds qdcoll((k-1)*d+1:k*d)];
    qdds = qdd(k);
    
    % loop over the collocation points
    for j=1:d
        % Expression of the state derivatives at the collocation points
        qp = qs_coll*C(:,j+1);
        qdp = qds_coll*C(:,j+1);         

        opti.subject_to(h*qds_coll(:,j+1) - qp == 0); % positions
        opti.subject_to(h*qdds  - qdp == 0);
    end
    
    % Feedback from reflexes
    if k > 7
        qs_delayed = q(k-7);
        qds_delayed = qd(k-7);
    else
        qs_delayed = q(1);
        qds_delayed = qd(1);
    end
    Tr  = (0.5*tanh(uR*qs_delayed  +   udotR*qds_delayed)+0.5)*(uR*qs_delayed  +   udotR*qds_delayed);
    Trh = (0.5*tanh(uR*-qs_delayed  +   udotR*-qds_delayed)+0.5)*(uR*qs_delayed  +   udotR*qds_delayed);
    
    % Limit torque
    theta_ref = (knee_r_range(2) - knee_r_range(1))/2;
    qrel = qs - (knee_r_range(1) + knee_r_range(2))/2;
    Tlim = -exp(klim*(qrel-theta_ref)) + exp(klim*(-qrel-theta_ref));
    
    % Skeleton dynamics
    opti.subject_to(qdds == (-uM*g*uL*cos(qs)- uDamp*qds + uTb + Tlim + Tr + Trh )./I);
    
    if k>= ki_StartObj
        qj = (qs-q_exp_interp(k)).^2;
        J = J + qj;
    end
    
    % state continuity at mesh transition
    opti.subject_to(q(k+1) - qs_coll*D == 0);
    opti.subject_to(qd(k+1)- qds_coll*D == 0);
end

opti.minimize(J);

opti.solver('ipopt',options);
sol = opti.solve();

%% Solutions
qsol      = sol.value(q)';
qdsol     = sol.value(qd)';
uDampsol  = sol.value(uDamp)
uTbsol    = sol.value(uTb)
uMsol     = sol.value(uM)
uLsol     = sol.value(uL)
uRsol     = sol.value(uR)
udotRsol  = sol.value(udotR)
uRhsol    = sol.value(uRh)
udotRhsol = sol.value(udotRh)

acc       = sum(abs(qsol(1:end-1)- q_exp_interp)) %Accuracy between simulated and experimental data 

%% Plot figure

figure
plot(tvect,qsol(1:end-1))
hold on
plot(tvect,q_exp_interp,'--k')
hold on
legend({'Casadi','Experimental'})

