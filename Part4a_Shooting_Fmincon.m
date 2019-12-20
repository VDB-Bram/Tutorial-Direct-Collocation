%% Optimization using shooting and Fmincon 

% This script explains how you can optimize the daming and baselinetorque
% to track recorded data of the pendulum test.

clear all; 

%% add functions to the matlab path
addpath(fullfile(pwd,'functions'));


%% Load pendulum properties (see Part1_Input.m for more information)

% Input
data  = load('DataPendulum.mat');
q_exp = data.data(:,2)*pi/180;
t_exp = data.data(:,1);
t_span = [t_exp(1) t_exp(end)]';

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

%% Optimization
    
% Optimization of Tb,B to minimize the difference between predicted
% and experimental values 

% initial guess of the optimization parameters
ig = [0 0];               % Tb B 

% upper and lower bounds on the optimization parameters
lb = [-10 -10] ;     % upper bound Tb B 
ub = [10  10] ;      % lowe rbound Tb B  

% define integration scheme for the forward simulation
dt      = 0.005;                    % time step
tvect   = t_span(1):dt:t_span(2);   % time vector
N       = length(tvect);

%  interpolate experimenatal data at discr. time using spline
qspline = spline(t_exp,q_exp);  % spline fit
[q_exp,qdot_exp] = SplineEval_ppuval(qspline,tvect,1); % get angles and velocities

% Define the objective function (myobj_shooting).
% Open the myobj_shooting to function to gain insight in how the shooting
% optimization is implemented.
F_handle  = @(z)myobj_shooting(z,tvect,dt,x0,q_exp,params);
[z, fval] = fmincon(F_handle,ig,[],[],[],[],lb,ub);

% Get the optimized variables
Tb = z(1);
B  = z(2);

% print results to screen
disp('Results:');
disp(['Optimal base line torque: ' num2str(Tb) 'Nm']);
disp(['Optimal damping coefficient ' num2str(B) 'Nms']);


%% Check simulation results using a forward simulation
   
params.Tb = Tb;
params.B  = B;


% Forward integration using time marching
dt      = 0.005;                    % time step
tvect   = t_span(1):dt:t_span(2);  % time vector
n       = length(tvect);           % number of points
x       = zeros(2,n);              % states
x(:,1)  = x0;                  % set initial state
for i=1:n-1
    % evaluate pendulum dynamics (using the function StateDerivative, note
    % that you already used this funciton in part 2C)
    xd = StateDerivative(tvect(i),x(:,i), params);
    % integration step (backward euler)
    x(:,i+1) = xd.*dt+x(:,i);
end

% get output
qsim      = x(1,:)';    % simulated angles
qdsim     = x(2,:)';    % simulated velocities
tsim      = tvect;
   
% plot output
figure();
plot(tvect,qsim); hold on;
plot(tvect,q_exp,'--k');
legend({'Q: Shooting','Q: Experimental'})
title('Shooting with optimized Tb and B ')
xlabel('Time [s]'); 
ylabel('Angle [rad]');


% save the results in the folder Results.
% Note that there is a script in the folder Results (CompareResults.m) that
% you can use at the end of the tutorial to compare shooting approach with
% the collocation approach.
save(fullfile(pwd,'Results','Shooting.mat'),'Tb','B','qsim','qdsim','tsim','q_exp');

