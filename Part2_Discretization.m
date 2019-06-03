%% Part 2 : Discretization - Backward Euler %%
clear all
close all
clc

% Input
data  = load('DataPendulum.mat');
q_exp = data.data(:,2)*pi/180;
t_exp = data.data(:,1);

m  = 2.3351;
lc = 0.2367;
gv = 9.81; 
RG = lc*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc;

% Torques
Tb = 2;
B  = 0.2;              % B will be used to calculate Td(damping) Td=-B*qdot

%% Part 2A : Time marching (for loop)

% Initial state
q0     = q_exp(1);
qdot0  = 0;
q_init = [q0 qdot0];
    
% Time info
dt       = 0.0001;                  % time step for backward euler
t_end    = 10;                      % final simulation time
t_vect   = 0:dt:t_end;              % vector with all the time steps
qt_store = nan(length(t_vect),2);   % pre-allocate matrix with predicted states (computation speed)

for i=1:length(t_vect)
    
    if i==1 % use initial state 
        q     = q0;
        qdot  = qdot0;
    else    
        q     = qi;    
        qdot  = qdoti;  
    end
    
    T         = Tb + B*qdot;
    qddot     = (-m*gv*lc*cos(q) - T)./I  ;
    
    qi        = qdot*dt+q;
    qdoti     = qddot*dt +qdot;
    
    % store the predicted state
    qt_store(i,1) = qi;
    qt_store(i,2) = qdoti;
    
end

% plot the predicted angle and angular velocity
figure();
subplot(1,2,1);
plot(t_vect,qt_store(:,1));xlabel('Time [s]'); ylabel('Angle [rad]');
hold on 
legend({'Q: simulated'})
title('Forward simulation of knee angle')
subplot(1,2,2);
plot(t_vect,qt_store(:,2)); xlabel('Time [s]'); ylabel('Angular velocity [rad/s]');
hold on 
legend({'Qdot: simulated'})
title('Forward simulation of knee angular velocity')

%% Part 2B: Implicit (fsolve)
% Initial state
q0     = q_exp(1);
qdot0  = 0;
q_init = [q0 qdot0];

params.m  = m;
params.lc = lc;
params.g  = 9.81;
params.I  = I;
params.Tb = Tb;
params.B  = B;

dt       = 0.01;                  % time step for backward euler
t_end    = 10;                      % final simulation time
t_vect   = 0:dt:t_end;  

x0 = zeros(length(t_vect),2);

F_handle = @(x)mycon_fsolve(x,dt,q_init,params)
x = fsolve(F_handle,x0);
q = x(:,1);
qdot = x(:,2);

figure()
plot(t_vect,q)

%% part 2C: Ode 

% Initial state
q0    = q_exp(1);
qdot0 = 0;
q_init = [q0 qdot0]';

% Time info
dt       = 0.0001;                  % time step for backward euler
t_end    = 10;                      % final simulation time
t_vect   = 0:dt:t_end;              % vector with all the time steps
t_span   = t_vect;

params.m  = m;
params.g  = gv;
params.lc = lc;
params.I  = I; 
params.Tb = Tb;
params.B  = B;
 
options = odeset('InitialStep',0.01,'MaxStep',0.01);

% Ode
[tM,qM] = ode23(@qdotfunctie_discretization, t_span, q_init,options, params);

figure()
plot(tM,qM(:,1))
xlabel('Time [s]'); 
ylabel('Angle [rad]');
legend({'Q: simulated'})
