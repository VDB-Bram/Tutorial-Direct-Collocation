
%% Load pendulum properties (see Part1_Input.m for more information)

% Input
data  = load('DataPendulum.mat');
q_exp = data.data(:,2)*pi/180;
t_exp = data.data(:,1);
t_span = [t_exp(1) t_exp(end)]';

% pendulum properties
m  = 2.3351;            % mass of the pendulum
lc = 0.2367;            % distance knee - COM tibia
g  = 9.81;              % gravitational constant
RG = lc*0.416;          % Radius of gyration (Winter 2009)
I  = m*RG*RG + m*lc*lc; % Inertia of lower limb +foot

% Initial state of the pendulum
x0 = [q_exp(1); 0];  % [joint angle initial frame  -  zero velocity ]
    

%% discretised time

% get the time vector for the simulation
dt       = 0.005;
tvect     = t_span(1):dt:t_span(end); 
N        = length(tvect);

% interpolate experimenatal data at discr. time using spline
qspline = spline(t_exp,q_exp);  % spline fit
[q_exp,qdot_exp] = SplineEval_ppuval(qspline,tvect,1); % get angles and velocities


%% formulate OCP

% import casadi libraries
import casadi.*;

% Initialise opti structure
opti = casadi.Opti();

% Discretized states
q 	= opti.variable(1,N);
qd 	= opti.variable(1,N);

% static opt parameters
Tb  = opti.variable(1);
B   = opti.variable(1);

% bounds  
opti.subject_to(-10 < Tb < 10);			% bound on baseline torque
opti.subject_to(-10 < B < 10);			% bound on damping coefficient
opti.subject_to(-4*pi < q <4*pi);	% bound on joint angles
opti.subject_to(-300 < qd < 300);		% bound on angular velocites

% bound on initial state
opti.subject_to(q(1)  == x0(1));
opti.subject_to(qd(1) == x0(2));

% guess
opti.set_initial(q, q_exp);
opti.set_initial(qd,qdot_exp);
opti.set_initial(Tb, 0);
opti.set_initial(B, 0);


%--------------------------
% using vector formulation
%--------------------------

% IP dynamics
qdd = (-m*g*lc*cos(q))/I - Tb/I -B*qd/I ; 	

% backward euler
opti.subject_to(qd(1:N-1)*dt +q(1:N-1) == q(2:N));
opti.subject_to(qdd(1:N-1)*dt +qd(1:N-1) == qd(2:N));

% -------------------------------------
% using for loop (easier to understand)
%--------------------------------------
% 
% for i = 1:N-1
% 	% select state
% 	qs = q(i);
% 	qds = qd(i);
% 	
% 	% dynamics
% 	qdd = (-m*g*lc*cos(q))/I - Tb/I -B*qd/I;
% 
% 	% state continuity at mesh transition (i.e. backward euler)
%     opti.subject_to(qds*dt+qs == q(i+1));
%     opti.subject_to(qdd*dt+qds == qd(i+1));
% 	
% end


% objective function
qerror = q - q_exp;
J      = sumsqr(qerror);
opti.minimize(J);

% options for IPOPT
if exist('options','var')
    clear options
end
options.ipopt.tol = 1*10^(-6);
options.ipopt.linear_solver = 'mumps';
options.ipopt.hessian_approximation = 'limited-memory';

% Solve the OCP
opti.solver('ipopt',options);
sol = opti.solve();

%% plot results
qsim  = sol.value(q)';
qdsim = sol.value(qd)';
B     = sol.value(B);
Tb    = sol.value(Tb);
tsim  = tvect;

figure()
plot(tsim,qsim,'--r'); hold on;
plot(tsim,q_exp, '--k')
legend({'Q: DC with Fmincon','Q: Experimental'})
xlabel('Time [s]'); 
ylabel('Angle [rad]');
title('DC - opti: optimized Tb and B ');

save(fullfile(pwd,'Results','Collocation_opti.mat'),'qsim','qdsim','B','Tb','tsim','q_exp');


