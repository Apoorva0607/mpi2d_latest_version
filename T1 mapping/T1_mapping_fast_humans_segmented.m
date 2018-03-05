clear all
% close all
clc




load /v/raid1a/apedgaon/MRIdata/Cardiac/Prisma/P120616/Processing_t1_seg/Pre_Interp_seg_3D_bin2_SRT1_MID_88.mat
 a=abs(double(new_reduced_non_k_space));
% a=double(sos_imgs);
load /v/raid1a/apedgaon/MRIdata/Cardiac/Prisma/P120616/Processing_t1_seg/Pre_Interp_seg_3D_bin2_SRT2_MID_88.mat
 b=abs(double(new_reduced_non_k_space));
 
 load /v/raid1a/apedgaon/MRIdata/Cardiac/Prisma/P120616/Processing_t1_seg/Pre_Interp_seg_3D_bin2_SRT3_MID_88.mat
 c=abs(double(new_reduced_non_k_space));
% b=double(sos_imgs);


inpath='/v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P111716/DicomData/';
instring='cine_tf2d12*';
outpath=['/v/raid1a/apedgaon/MRIdata/Cardiac/Prisma/P111716/Processing_T1/Dicoms/Pre LV/'];
seriesNum=12345;
seriesDesc='3D T1 maps LV_pre';
SRT1=200;
SRT2=500;


for i=1:30
   i 
%    img_series=double(abs(squeeze(new_reduced_non_k_space(:,:,i,:))));
%    S1=medfilt2(mean(img_series(:,:,2:2:end),3),[5 5]);
%    S2=medfilt2(mean(img_series(:,:,3:2:end),3),[5 5]);
%    A=medfilt2(img_series(:,:,1),[5 5]);
   

%  
%     S1=(a(:,:,i)).*1e7.*(0.2/0.4);
%     S2=(b(:,:,i)).*1e7.*(0.2/0.4);
%     A=c(:,:,i).*1e7;
load('bins_MID88.mat');

             S1=((a(:,:,i)).*1e7./0.35);
        S_1=((a(:,:,i)).*1e7./0.35)/bin2_srt1(i);
    S2=(b(:,:,i)).*1e7./0.35;
    S_2=(b(:,:,i)).*1e7./0.35/bin2_srt2(i);
    
A=c(:,:,i).*1e7./0.3;
A1=c(:,:,i).*1e7./0.3/bin2_srt3(i);

%     S1=medfilt2(a(:,:,i),[5 5]);
%     S2=medfilt2(b(:,:,i),[5 5]);
%     A=medfilt2(pd_img(:,:,i),[5 5]);



   
   S1_A=1-(S_1./A1);
   S2_A=1-(S_2./A1);

   
   temp_S1=find(S1_A<0);
   S1_A(temp_S1)=0;
   temp_S2=find(S2_A<0);
   S2_A(temp_S2)=0;

   
   
   T1_temp_num=SRT2*SRT1;
   T1_temp_den=log(S1_A).*log(S2_A);
   
  temp_den=find(T1_temp_den<0);
  T1_temp_den(temp_den)=1; 
   
   T1_bin2_88(:,:,i)=sqrt((T1_temp_num./T1_temp_den));
   
    
end

%   T1map_dicoms(T1,inpath,instring,outpath,seriesNum,seriesDesc)