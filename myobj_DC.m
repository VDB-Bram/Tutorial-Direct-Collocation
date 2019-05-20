function [f_out] = myobj_DC(z,params,t)

%MYOBJ is the function dat you want to minimize 

%Input 
q       = z(3:2+params.N);
q_exp   = params.q_exp;

  
% Output 
f_out = (q - q_exp').^2;
f_out = sum(f_out);

end

