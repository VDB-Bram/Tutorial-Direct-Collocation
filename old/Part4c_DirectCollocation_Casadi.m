%% Direct Collocation - Casadi 
clear all
clc
close all
    
% Input
data      = load('DataPendulum.mat');
t_exp     = data.data(:,1);
q_exp     = data.data(:,2)*pi/180;
q_exp_dot = diff(q_exp)/0.01; 
q_exp_dot(end+1) = q_exp_dot(end); 
   
m   = 2.3351;
lc  = 0.2367;
gv  = 9.81; 
RG  = lc*0.416;          % Radius of gyration (Winter 2009)
I   = m*RG*RG + m*lc*lc;

% Time 
dt     = 0.01;
t_vect = 0:dt:(t_exp(end)- t_exp(1));  
% t_vect = t_exp(1):dt:t_exp(end)+dt; 
N      = length(t_vect);  

%% Load Casadi and linear solver
addpath('C:\software\linear_solver');   
Casadi = 'C:\software\Casadi'; addpath(genpath(Casadi));
import casadi.*;

%% Casadi

% setup polynomials
d = 3;                  % degree of polynomials
tau_root = [0 collocation_points(d, 'radau')];
C = zeros(d+1,d+1);     % Coefficients of the collocation equation
D = zeros(d+1, 1);      % Coefficients of the continuity equation
E = zeros(d+1, 1);      % Coefficients of the quadrature function

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
  E(j) = polyval(pint, 1.0);
end

% Start with an empty NLP
w   = {};
w0  = [];
lbw = [];
ubw = [];
J   = 0;
g   = {};
lbg = [];
ubg = [];

% model variables
q    = SX.sym('q');
qd   = SX.sym('qd');
x    = [q; qd];
Tb   = SX.sym('uTb');
B    = SX.sym('uB');

% Skeleton dynamics
xdot = [ qd; (-m*gv*lc*cos(q)- B*qd + Tb )./I]; 
f = Function('f', {x,B, Tb}, {xdot});             % continuous time dynamics

% Damping
uB = MX.sym('uB', 1);
w = {w{:}, uB};
lbw = [lbw; 0.01];
ubw = [ubw; 10];
w0 = [w0; 0.01]; 

% Tb
uTb = MX.sym('uTb', 1);
w = {w{:}, uTb};
lbw = [lbw; 0];
ubw = [ubw; 10];
w0 = [w0; 0]; 

% Initial conditions q and qdot 
Xk = MX.sym('X0', 2);
w = {w{:}, Xk};
lbw = [lbw; q_exp(1); 0];
ubw = [ubw; q_exp(1); 0];
w0 = [w0; q_exp(1); 0]; 


% Formulate the NLP
for k=0:N-1
   
    % State at collocation points
    Xkj = {};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], 2);
        w = {w{:}, Xkj{j}};
        lbw = [lbw; -300; -300];
        ubw = [ubw;  300;  300];
        w0 = [w0; 0; 0];            % guess set at zero
    end

    % Loop over collocation points
    Xk_end = D(1)*Xk;
    for j=1:d
       % Expression for the state derivative at the collocation point
       xp = C(1,j+1)*Xk;
       for r=1:d
           xp = xp + C(r+1,j+1)*Xkj{r};
       end

       % Append collocation equations
       [fj] = f(Xkj{j},uTb,uB);
       
       g    = {g{:}, dt*fj - xp};
       lbg  = [lbg; 0; 0];       % violation on dynamics should be zero
       ubg  = [ubg; 0; 0];       % violation on dynamics should be zero

       % Add contribution to the end state
       Xk_end = Xk_end + D(j+1)*Xkj{j};

       % Add contribution to quadrature function
       qj = (Xkj{j}(1)-q_exp(k+1)).^2;     
       J  = J + E(j+1)*qj*dt;
    end

       % New NLP variable for state at end of interval
       Xk = MX.sym(['X_' num2str(k+1)], 2);
       w = {w{:}, Xk};
       lbw = [lbw; -300;  -300];
       ubw = [ubw;  300;  300];
       w0 = [w0; q_exp(k+1); q_exp_dot(k+1)];
  
       % Add equality constraint
       g = {g{:}, Xk_end-Xk};
       lbg = [lbg; 0; 0];
       ubg = [ubg; 0; 0];
       
end

% Create an NLP solver
prob   = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
options.ipopt.linear_solver = 'ma57';
solver = nlpsol('solver', 'ipopt', prob,options);

% Solve the NLP
tic
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
toc
w_opt = full(sol.x);

% get the output
q_opt = w_opt(3:2*d:end);
qd_opt = w_opt(4:2*d:end);


%% Plot the output
figure()
plot(t_vect,q_opt(1:N)) 
hold on
plot(t_vect,q_exp,'--k')
hold on
legend({'Q: Simulated','Q: Experimental'})
xlabel('Time [s]'); 
ylabel('Angle [rad]');
title('DC - Casadi: optimized Tb and B ')
