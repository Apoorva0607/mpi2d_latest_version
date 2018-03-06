clear all
close all
clc
warning off

 dataPath='/v/raid1a/dlikhite/MRIdata/Cardiac/Prisma/P082616/ReconData/';
oldoptions=optimset('fmincon');
options=optimset(oldoptions,'TolFun',0.00000001,'Display','off') ;
inp_ser=[89,90,91];
% inp_ser=[145,144,143,142,140];
SRT=[400,600];
slice=6;
% load mask.mat
mask=ones(512,512);
ref_img=ones(512,512,3);
% try
% load('t1_series_moco.mat');    
%     
% catch
for i=3
    
file=([dataPath,'Pre_Interp_3D_MID_',int2str(inp_ser(i)),'.mat']);
load(file)

tmp_img=squeeze(double(1e10.*new_reduced_non_k_space(:,:,slice,:)));
tmp_img1=abs((tmp_img(:,:,1:15)));

t1_series(:,:,:,i)=tmp_img1;
ref_img(:,:,i)=t1_series(:,:,3).*ref_img(:,:,i);
end



%  load t1_series_moco.mat
% t1_series_tmp=t1_series(:,:,2:end,:);
% t1_series_tmp1=squeeze(mean(t1_series_tmp,3));
% for i=1:3
%     ref_img(:,:,i)=ref_img(:,:,i).*t1_series_tmp1(:,:,3);
% end
% 
% 
% 
% 
% t1_series_moco=ANTS_DMBR(t1_series_tmp1,ref_img,mask,0.1,0.5,0.5,1);      
% 
% t1_series_moco(:,:,3)=t1_series_tmp1(:,:,3);
% save('t1_series_moco.mat','t1_series_moco');

% 
% t1_series_moco_tmp=t1_series_moco(:,:,:);
% t1_map=zeros(512,512);
% 
% 
% A=[];
% B=[];
% nonlcon=[];
% x0(1)=double(1000.0);
% % x0(2)=0;
% LB(1)=double(0.0);
% UB(1)=double(5000.0);
% % LB(2)=double(0.0);
% % UB(2)=double(1e5);
% sin_a=sin(2*pi/180);
% sin_alpha=sin(12*pi/180);
% 
% 
% for i=273
%     i
%    for j=237
% 
%  
%  
%  A_scalar=t1_series_moco_tmp(i,j,3).*(sin_alpha/sin_a);
% curve=(squeeze(t1_series_moco_tmp(i,j,1:2)));
% 
% t1_map(i,j)= fmincon('T1_3param',x0,A,B,A,B,LB,UB,nonlcon,options,curve,A_scalar,SRT);
% 
% 
% % A_scalar(i,j)=x(2);
%     end
% end