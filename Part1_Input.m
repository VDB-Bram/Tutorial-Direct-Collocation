%% Part 1: Equations of motion %%
clear all
close all
clc

% Pendulum movement

%% Input values

m  = 2.3351;
lc = 0.2367;
gv = 9.81; 
RG = lc*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc;

% Torques
Tb = 2;
B  = 0.2;  % B will be used to calculate Td(damping) Td=B*qdot


%% Experimental data 

data  = load('DataPendulum.mat');
q_exp = data.data(:,2)*pi/180;
t_exp = data.data(:,1);

% plot

figure()
plot(t_exp,q_exp,'--k')
legend({'Q: Experimental'})
xlabel('Time [s]'); 
ylabel('Angle [rad]');
title ('Kinematics of pendulum movement (experimental)')
