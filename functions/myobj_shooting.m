function [f_out] = myobj_shooting(z,tvect,dt,x0,q_exp,params)

%MYOBJ evaluates the objective function of a tracking simulations in the
%pendulum test using a shooting approach
% output:
%   f_out   = Difference between q_exp (experimental data) and xM (predicted
%           values)

% intput:
%   z = vector with optimization variables
%       z(1)    = Tb
%       z(2)    = B
%   tvect   = discretised time
%   dt      = time step used for discritization
%   x0      = initial state (joint angles and velocities in radians)
%   q_exp   = experimental data at discretised time
%   params  = optimization parameters

% store the optimization variables (Tb and B) in the params structure
params.Tb = z(1);
params.B  = z(2);

% Forward integration using time marching
n       = length(tvect);           % number of points
x       = zeros(2,n);              % states
x(:,1)  = x0;                  % set initial state
for i=1:n-1
    % evaluate pendulum dynamics (using the function StateDerivative, note
    % that you already used this funciton in part 2C)
    xd = StateDerivative(tvect(i),x(:,i), params);
    % backward euler
    x(:,i+1) = xd.*dt+x(:,i);
end

% evaluate difference between measured and simulated kinematics
qError = x(1,:) - q_exp;
f_out = sumsqr(qError);

end

