function [Fsum]=T1_3param(x,curve,A_scalar,SRT)

warning off

T1=x(1);
% A_scalar=x(2);
errsum = 0.0;

for i=1:length(SRT)
curve_est(i)=A_scalar.*(1-exp(-SRT(i)/T1));        
    
end

% plot(curve)
% hold on
% plot(curve_est)
% pause(0.1)

err = (curve' - curve_est);
      
errsum = errsum + sum(err.*err);

Fsum=errsum;