function gd_curves=mpi_si2gdconc_test1(Mz0,curves,slice_sat_no,sin_a,sin_asmall,TD,nor,flag_PD_end)
TD=25;
nor=24;


  curves=[Mz0 curves(:,26:end)];
    curves1=curves(:,21:end);


 
  

for i=1:size(curves1,1)
   si0=mean(curves1(i,1:3));
       
   for j=1:size(curves1,2)
       
      curves2(i,j)=(curves1(i,j)-si0)./si0; 
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag_multihance=1;

Rbld=zeros(1,size(curves2,1));  %prohance=3.8;  % from donahue review '97   %multihance=6.3 .....Issue: Volume 41(3), March 2006, pp 213-221
if flag_multihance==1
    Rbld(1)=6.3;
    Rbld(2:end)=3.8;
else
    Rbld(1:end)=3.8;
end

bld_T1_0=1.666; % larsson 2001    %1.379;
 % larsson 2001     %1.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Trec=(TD+(slice_sat_no-1).*2.24*nor).*0.001;


for i=1:1:size(curves2,1)



for l=1:1:size(curves2,2)
 num=-Trec;
 den=log((curves2(i,l)+1).*exp(-Trec./bld_T1_0)-curves2(i,l));
 t1_unscaled(l)=num./den;
 
 
end

for k=1:1:length(t1_unscaled)
    C_unscaled(k)=((1/t1_unscaled(k))-(1/bld_T1_0))/Rbld(i);
end

 scale=sum(C_unscaled(11:12)./2);
 gd_curves(i,:)=(C_unscaled);

end
 %gd_curves(find(gd_curves<0))=0;
%  gd_curves=[zeros(size(gd_curves,1),20),gd_curves];
plot(gd_curves')
hold on

return;