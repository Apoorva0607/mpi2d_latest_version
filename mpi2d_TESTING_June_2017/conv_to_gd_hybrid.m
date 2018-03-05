function Fsum=conv_to_gd_hybrid(x,Mz,Trec,value,n,Tr,alpha)



t1=x(1);

% sin_a=sin(x(3)*pi/180);

Tr=Tr*0.001; 

sin_a=sin(alpha.*pi/180);
cos_a=cos(alpha.*pi/180);


a1=cos_a*exp(-Tr/t1);

b1=1-exp(-Tr/t1);

% n=ones(1,24);

% if count==1;
%  for i=1:length(n)
% S_tmp(i)=(Mz.*sin_a)*((1-exp(-Trec/1.6))*(a1^(n(i)-1))+b1*((1-a1^(n(i)-1))/(1-a1)));
% 
% end
% S=mean(S_tmp);   
% 
% 
% scale=value./S;
% save('scale.mat','scale');
% end
% 
% load('scale.mat')


for i=1:length(n)
S_tmp(i)=(Mz.*sin_a)*((1-exp(-Trec/t1))*(a1^(n(i)-1))+b1*((1-a1^(n(i)-1))/(1-a1)));

end

% S_tmp=(1.94.*Mz0*sin_a./sin_a_small)*(1-exp(-Trec/t1));

S=mean(S_tmp);
Fsum=abs((S)-value);


