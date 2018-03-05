function [t1_value,t1_value1,C_unscaled,C_unscaled1] = test_gd_hybrid(M,Trec,S,n,X0,bld_T1_0,Rbld,scale_pd,Tr,flip_pd,alpha)


warning off all
oldoptions = optimset('Display', 'off') ;

options=optimset(oldoptions,'TolFun',0.0001,'TolX',0.0001,'MaxIter',1e4) ;
A=[];
B=[];
nonlcon=[];

LB(1)=-10;
UB(1)=10;

% if scale_pd==0
% LB(2)=0;
% UB(2)=10000;
% clear X0
% X0=5000;
% count=1;    
% for i=1:5
% [x(count)]=fminsearch('conv_to_gd_hybrid_Mz',X0,options,Trec,(S(i)),n);
% %[x(count)]=fmincon('conv_to_gd_hybrid_Mz',X0,A,B,A,B,LB,UB,nonlcon,options,Trec,(S(count)),n);
% count=count+1;
% X0=mean(x);
% end
% end





%M=max(x);
M=(M/sin(flip_pd*pi/180));


for count=1:length(S)


% LB(3)=9;
% UB(3)=11;

value=S(count);

    [x]=fmincon('conv_to_gd_hybrid',X0,A,B,A,B,LB,UB,nonlcon,options,M,Trec,value,n,Tr,alpha);
%     [x]=fminsearch('conv_to_gd_hybrid',X0,options,M,Trec,(value),n,Tr,alpha);
X0=x;
t1_value(count)=x(1);
t1_value1(count)=0;
% Mz(count)=x(2);

C_unscaled(count)=((1/t1_value(count))-(1/bld_T1_0))/Rbld;
C_unscaled1(count)=((1/t1_value1(count))-(1/bld_T1_0))/Rbld;
end