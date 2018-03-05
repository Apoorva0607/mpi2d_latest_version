function Fsum=t1_fun(x,S1,S2,S_PD,SRT1,SRT2)



T1=x(1);


%   Fsum=abs((S1/S2)-((1-exp(-SRT1/T1))/(1-exp(-SRT2/T1))));

  Fsum=abs(((S1*sin(2*pi/180))/(S_PD*sin(12*pi/180)))-((1-exp(-SRT1/T1))));
%   
%   Fsum2=abs(((S2*sin(2*pi/180))/(S_PD*sin(12*pi/180)))-((1-exp(-SRT2/T1))));
%   
%   Fsum=abs(Fsum1-Fsum2);