clear all
close all
clc
warning off all
oldoptions = optimset('Display', 'off') ;

options=optimset(oldoptions,'TolFun',1e-30,'TolX',1e-30,'MaxIter',1e30) ;
A=[];
B=[];
nonlcon=[];
LB=0;
UB=4000;
SRT1=400;
SRT2=600;
x0=700;

load t1_series_moco.mat

t1_series_moco_tmp=double(t1_series_moco(150:350,150:350,:));


parfor i=1:200
    i
    for j=1:200
        
        
        S1=1.43.*t1_series_moco_tmp(i,j,1);
        S2=1.43.*t1_series_moco_tmp(i,j,2);
        S_PD=t1_series_moco_tmp(i,j,3);
        
        if S1>0 && S2>0
            t1_map(i,j)=fmincon('t1_fun',x0,A,B,A,B,LB,UB,nonlcon,options,S1,S2,S_PD,SRT1,SRT2);
        else
            t1_map(i,j)=0;
        end
        
        
        
    end
end