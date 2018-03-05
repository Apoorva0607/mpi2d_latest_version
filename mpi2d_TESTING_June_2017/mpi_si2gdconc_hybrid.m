function [gd_curves T1]=mpi_si2gdconc_hybrid(Mz,curves_tmp,slice_sat_no,sliceNum,seriesNum)

scale_pd=1;

TD=20;  
Tr=1.6;
flip_pd=2;
alpha=10;
bld_T1_0=1.6; % larsson 2001    %1.379;
X0(1)=bld_T1_0;


% a=round(seriesNum/100);
% a1=a/38;
% nor=(a1-1)*36+1:a1*36; 
 %nor=241:276;
%   nor=89:124;
% nor=241:276;
% nor=106:165;
nor=1:23;

numPreContrast=3;
numSkip=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if scale_pd==1

% if sliceNum==2
%     load('Output/curves.study4301.slice1.mat')
% elseif sliceNum==4
%     load('Output/curves.study4401.slice2.mat');
% else
%     load('Output/curves.study4501.slice3.mat');
% end

for i=1:size(curves_tmp,1)
% M(i)=max(Mz(i,2:end-2));
M(i)=Mz(i);
end
clear curves



curves=curves_tmp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    X0(2)=5000;
    M=zeros(1,7);
    curves1=[curves_tmp];
    
    for i=1:1:size(curves1,1)
        
    init_tiss(i)=sum(curves1(i,1+numSkip:numPreContrast+numSkip))/numPreContrast;
          
    end
    
for i=2:1:size(curves1,1)
curves(i,:)=(curves1(i,:)./((init_tiss(i)))).*mean(mean(init_tiss(2:end)));
end
curves(1,:)=curves1(1,:);

% scale_local=mean(mean(curves(2:end,5:8)));
% load('Output/scale_glo.mat');
% 
%  curves=(curves./scale_local).*scale_glo;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag_multihance=0;

Rbld=zeros(1,size(curves,1));  %prohance=3.8;  % from donahue review '97   %multihance=6.3 .....Issue: Volume 41(3), March 2006, pp 213-221
if flag_multihance==1
    Rbld(1)=6.3;
    Rbld(2:end)=3.8;
else
    Rbld(1:end)=3.8;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% if slice_sat_no==1;
% Trec=(TD).*0.001;
% else
%  Trec=(TD+2.2*nor(1)).*0.001;  
% end
   Trec=(TD+(slice_sat_no-1).*Tr*nor(end)).*0.001;

 
for i=1:1:size(curves,1)

disp('Wait')

[t1_value_a,t1_value_b,C_unscaled,C_unscaled1]=test_gd_hybrid(M(i),Trec,curves(i,:),nor,X0,bld_T1_0,Rbld(i),scale_pd,Tr,flip_pd,alpha);
 
% C_unscaled=real(C_unscaled);
%    C_unscaled(find(C_unscaled<0))=0;
%  C_unscaled(find(C_unscaled>10))=0;

T1(i,:)=t1_value_a;
scale=sum(C_unscaled(3:6)./4);

gd_curves(i,:)=(C_unscaled)-scale;


end


return;
