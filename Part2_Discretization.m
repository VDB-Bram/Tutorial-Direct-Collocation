%% Part 2 : Discretization - Backward Euler %%

%% Part 2A : Time marching (for loop)

% Initial state
q0    = q_exp(1);
qdot0 = 0;
    % or in matrix (x0=[q_exp(1) 0])
    
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
    qddot     = (-m*g*lc*cos(q) - T)./I  ;
    
    qi        = qdot*dt+q;
    qdoti     = qddot*dt +qdot;
    
    % store the predicted state
    qt_store(i,1) = qi;
    qt_store(i,2) = qdoti;
    
end

% plot the predicted angle and angular velocity
figure();
subplot(1,2,1);
plot(t_vect,qt_store(:,1));xlabel('Time [s]'); ylabel('angle [rad]');
subplot(1,2,2);
plot(t_vect,qt_store(:,2)); xlabel('Time [s]'); ylabel('Angular velocity [rad/s]');

%% Part 2B: implicit (fsolve)

%% part 2C: ode 
