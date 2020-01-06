%% Part 1: Pendulum test

% This script contains information about the properties of the shank and
% the measured motion during the pendulum test.

%% pendulum properties
m  = 2.3351;            % mass of the pendulum
lc = 0.2367;            % distance knee - COM tibia
gv = 9.81;              % gravitational constant
RG = lc*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc; % Inertia of lower limb +foot

% Torques
Tb = 2;                 % Baseline torque of 2 Nm
B  = 0.2;               % Damping coefficient (Nm s)

%% Experimental data 

% Load matfile with experimental data of the pendulum movement
data  = load('DataPendulum.mat');
q_exp = data.data(:,2)*pi/180;
t_exp = data.data(:,1);


%% Plot figure with measured joint angle as a function of time
figure()
plot(t_exp,q_exp,'--k')
legend({'Q: Experimental'})
xlabel('Time [s]'); 
ylabel('Angle [rad]');
title ('Kinematics of pendulum movement (experimental)')
