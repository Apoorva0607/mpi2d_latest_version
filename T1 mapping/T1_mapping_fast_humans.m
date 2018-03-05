clear all
% close all
clc




load /v/raid11/gadluru/MRIdata/Cardiac/Prisma/P093016/ReconData/mat_files/3D/T1_mapping/set1/Pre_Interp_3D_MID_70.mat
 a=abs(double(new_reduced_non_k_space));
% a=double(sos_imgs);
load /v/raid11/gadluru/MRIdata/Cardiac/Prisma/P093016/ReconData/mat_files/3D/T1_mapping/set2/Pre_Interp_3D_MID_70.mat
 b=abs(double(new_reduced_non_k_space));
% b=double(sos_imgs);
load pd_img

inpath='/v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P093016/DicomData/';
instring='cine_tf2d12*';
outpath=['/v/raid1a/dlikhite/MRIdata/Cardiac/Prisma/P093016/Processing_T1/Dicoms/Pre LV/'];
seriesNum=12345;
seriesDesc='3D T1 maps LV_pre';
SRT1=200;
SRT2=650;


for i=1:8
   i 
%    img_series=double(abs(squeeze(new_reduced_non_k_space(:,:,i,:))));
%    S1=medfilt2(mean(img_series(:,:,2:2:end),3),[5 5]);
%    S2=medfilt2(mean(img_series(:,:,3:2:end),3),[5 5]);
%    A=medfilt2(img_series(:,:,1),[5 5]);
   

 
    S1=(mean(squeeze(a(:,:,i,2:end-1)),3));
    S2=(mean(squeeze(b(:,:,i,2:end-1)),3));
    A=abs(a(:,:,i));


%     S1=medfilt2(a(:,:,i),[5 5]);
%     S2=medfilt2(b(:,:,i),[5 5]);
%     A=medfilt2(pd_img(:,:,i),[5 5]);



   
   S1_A=1-(S1./A);
   S2_A=1-(S2./A);

   
   temp_S1=find(S1_A<0);
   S1_A(temp_S1)=0;
   temp_S2=find(S2_A<0);
   S2_A(temp_S2)=0;

   
   
   T1_temp_num=SRT2*SRT1;
   T1_temp_den=log(S1_A).*log(S2_A);
   
  temp_den=find(T1_temp_den<0);
  T1_temp_den(temp_den)=1; 
   
   T1(:,:,i)=sqrt((T1_temp_num./T1_temp_den));
   
    
end

  T1map_dicoms(T1,inpath,instring,outpath,seriesNum,seriesDesc)