function op_img=registration_standalone(in_img,mask,iseries,isub,slice,old_flag)

PD_keep=0;
pd_end=1;
rangex=size(in_img,1);
rangey=size(in_img,2);

if old_flag==0
    step=0.2;
    a1=0.5;
    a2=0.5;
else
    step=0.1;
    a1=1;
    a2=1;
end
    


if PD_keep==1
    cinemri1_tmp=in_img(:,:,pd_end+1:end);
    cinemri1_pd=in_img(:,:,1:pd_end);
else
    cinemri1_tmp=in_img(:,:,pd_end+1:end);
end




run_registration(cinemri1_tmp,rangey,rangex,25, 2,2,2,iseries+isub);

temp=strcat('Output/',int2str(iseries+isub),'/cinemri_curv_fit1.mat');
load(temp);
temp=strcat('Output/',int2str(iseries+isub),'/buffer_cinemri.mat');
load(temp);


    cinemri2=ANTS_DMBR(buffer_cinemri,cinemri_curv_fit,mask,step,a1,a2,old_flag);

% run_registration(cinemri2_tmp,rangey,rangex,25, 3,3,3,iseries+isub+10);
%     temp=strcat('Output/',int2str(iseries+isub+10),'/cinemri_curv_fit1.mat');
%   load(temp);
% 
%   temp=strcat('Output/',int2str(iseries+isub),'/buffer_cinemri.mat');
%  load(temp);
%  
%  
%  cinemri2=ANTS_DMBR(buffer_cinemri,cinemri_curv_fit,step,a1,a2);
 
 !rm -f *mhd
 !rm -f *raw
 !rm -f *gz
 !rm -f *txt
 
  if PD_keep==1 
 op_img=zeros(size(cinemri2,1),size(cinemri2,2),size(cinemri2,3)+pd_end);
  op_img(:,:,1:pd_end)=cinemri1_pd(:,:,:);
 op_img(:,:,pd_end+1:end)=cinemri2;
 else
 op_img=cinemri2;
  end
 
  
  
 
  
  