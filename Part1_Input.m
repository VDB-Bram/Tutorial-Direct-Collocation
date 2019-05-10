%% Part 1: Equations of motion %%

% Pendulum movement

%% Input values

m  = 2.318;
lc = 0.2368;
g  = 9.81; 
RG = lc*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc;

% Torques
Tb = 2;
B  = 0.2;  % B will be used to calculate Td(damping) Td=B*qdot


%% Experimental data 

data  = load('DataPendulum.mat');
q_exp = data.data(:,2);
t_exp = data.data(:,1);


