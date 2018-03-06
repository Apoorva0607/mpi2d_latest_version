clear all
clc
close all

oldoptions=optimset('fmincon');
options=optimset(oldoptions,'TolFun',0.00000001,'Display','off') ;
A=[];
B=[];
nonlcon=[];
x0(1)=double(700.0);
LB(1)=double(0);
UB(1)=double(2000.0);



load t1_series_41.mat
t1_series_tmp=double(t1_series(:,:,:));
SRT=[202,400];

for i=1:size(t1_series_tmp,1)
    i
    parfor j=1:size(t1_series_tmp,2)
        A_scalar=t1_series_tmp(i,j,1);
        curve=squeeze(t1_series_tmp(i,j,2:3));
        
        if A_scalar>1
            t1_map(i,j)= fmincon('T1_3param',x0,A,B,A,B,LB,UB,nonlcon,options,curve,A_scalar,SRT);
        else
            t1_map(i,j)=0;
        end
        
    end
end