function [f_out] = myobj_shooting(z,q_exp,t_span,q_init,t_exp,params)

%MYOBJ f_out is difference between predicted and experimental values 
%   f_out   = Difference between q_exp (experimental data) and xM (predicted
%           values)
%   z(1)    = Tb
%   z(2)    = B
%   q_exp   = experimental data
%   t_span  = start and end point of experimental data
%   q_init  = Initial angle to start (radians)
%   t_exp   = time of experimetnal data 

t_exp = t_exp';

% Optimization values 
%   options   = [];
    options   = odeset('InitialStep',0.01,'MaxStep',0.01);
    params.Tb = z(1);
    params.B  = z(2);
  
% Forward integration 
  [tM,qM] = ode23(@qdotfunctie_shooting, t_span, q_init, options, params); 

% Interpolation 
  qMinterp = interp1(tM,qM(:,1),t_exp);
   
% F_out
  f_out = (qMinterp - q_exp).^2;
  f_out = trapz(t_exp,f_out);

end

