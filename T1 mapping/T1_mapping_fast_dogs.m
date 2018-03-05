clear all
clc
close all
warning off
oldoptions=optimset('fmincon');
options=optimset(oldoptions,'TolFun',0.00000001,'Display','off') ;
A=[];
B=[];
nonlcon=[];
x0(1)=double(700.0);
LB(1)=double(0);
UB(1)=double(2000.0);


inpath=[pwd,'/Pre_PD_SA_40'];
instring='IM*';
outpath=[pwd,'/Dicom/'];
seriesNum=20001;
seriesDesc='3D T1 maps pre LA';


list_dir=dir('Pre*');

for dir_count=1:length(list_dir)

    cd(list_dir(dir_count).name);
    list_dicoms=dir('*.dcm');
    
    for dicom_count=1:length(list_dicoms)
        
        t1_series(:,:,dir_count,dicom_count)=double(dicomread(list_dicoms(dicom_count).name));
        
        
    end
    cd ..
end

SRT1=200;
SRT2=800;


for i=1:96
   
   img_series=double(squeeze(abs(t1_series(:,:,:,i))));
   S1=img_series(:,:,2);
   S2=img_series(:,:,3);
   A=img_series(:,:,1);
  
   S1_A=1-(S1./A);
   S2_A=1-(S2./A);

% S1_A=A-S1;
% S2_A=A-S2;
%    
%    temp_S1=find(S1_A<0);
%    S1_A(temp_S1)=0;
%    temp_S2=find(S2_A<0);
%    S2_A(temp_S2)=0;
   
   
   
   T1_temp_num=SRT2*SRT1;
   T1_temp_den=log(S1_A).*log(S2_A);
   
%   temp_den=find(T1_temp_den<0);
%   T1_temp_den(temp_den)=1; 
   
   T1(:,:,i)=(sqrt(T1_temp_num./T1_temp_den));
      
   
end
 
   
    T1map_dicoms(abs(T1),inpath,instring,outpath,seriesNum,seriesDesc)