function gd_curves=mpi_si2gdconc(Mz0,curves,slice_sat_no,sin_a,sin_asmall,TD,nor,flag_PD_end,pd_end)
TD=25;
sin_a=10;   %HARDCODED FOR REPEATIBILITY STUDY
nor=24;
sin_asmall=2.5;


if flag_PD_end==0
    curves=[Mz0 curves(:,pd_end+1:end)];
    Mz0=Mz0(:,1:pd_end);
    curves1=curves(:,pd_end+1:end);
    
else
    Mz0=Mz0(:,1:2);
    curves1=curves(:,1:end-2);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag_multihance=0;

Rbld=zeros(1,size(curves1,1));  %prohance=3.8;  % from donahue review '97   %multihance=6.3 .....Issue: Volume 41(3), March 2006, pp 213-221
if flag_multihance==1
    Rbld(1)=6.3;
    Rbld(2:end)=3.8;
else
    Rbld(1:end)=3.8;
end

bld_T1_0=1.666; % larsson 2001    %1.379;
 % larsson 2001     %1.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Trec=(TD+(slice_sat_no-1).*2.2*nor).*0.001;
sin_a=sin(sin_a.*pi/180);
sin_asmall=sin(sin_asmall.*pi./180);

for i=1:1:size(curves1,1)

 Mz0_1=max(Mz0(i,:))./(sin_asmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for l=1:1:size(curves1,2)
        num=-Trec;
        den=log(1-(curves1(i,l)./(Mz0_1*sin_a)));
        t1_unscaled(l)=num./den;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     for l=1:size(curves1,2)
%        eqn=Mz0_1.*sin_a((1-exp(-Trec/t1))*(cos_a.*exp(-Trec/t1)).^(n-1) + (1-exp(-Trec/t1)).*( (1-(cos_a.*exp(-Trec/t1)).^(n-1))/(1-(cos_a.*exp(-Trec/t1))))) == curves1(i,l); 
%         
%        t1_unscaled_test1(l)=solve(eqn,t1); 
% 
%     end
% 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    for k=1:1:length(t1_unscaled)
        C_unscaled(k)=((1/t1_unscaled(k))-(1/bld_T1_0))/Rbld(i);
        
    end

 scale=sum(C_unscaled(2:5)./4);

 gd_curves(i,:)=(C_unscaled)-scale;

end

 %gd_curves(find(gd_curves<0))=0;
 
%  gd_curves=[zeros(size(gd_curves,1),pd_end),gd_curves];
return;
