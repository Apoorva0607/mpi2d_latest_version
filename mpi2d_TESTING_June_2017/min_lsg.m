
function mr_sh=min_lsg(rangex,rangey,slice,MBR_iter_no)

 temp=strcat('Output/',int2str(slice),'/cinemri.mat');
 load(temp);
crop_cinemri=cinemri;


temp=strcat('Output/',int2str(slice),'/cinemri_curv_fit',int2str(MBR_iter_no),'.mat');
 load(temp);
crop_cinemri_curv_fit=cinemri_curv_fit;

tmpp=mean(mean(mean(cinemri)));
cinemri=cinemri-tmpp;

temp_var=zeros(size(crop_cinemri_curv_fit));

mr_sh=zeros(size(cinemri_curv_fit,3),2);

clear m2 sh;
for i=1:size(cinemri,3)    
    [m2 sh]=mpi_leastsquares(crop_cinemri,crop_cinemri_curv_fit(:,:,i),crop_cinemri(:,:,i),i,rangex,rangey,slice);        
    temp_var(:,:,i)=m2;    
    mr_sh(i,:)=sh;    
end

 temp=strcat('Output/',int2str(slice),'/mr_sh.mat');
save(temp,'mr_sh');

 temp=strcat('Output/',int2str(slice),'/buffer_cinemri.mat');
load(temp);

tmp_imgs=zeros(size(buffer_cinemri));
for i=1:size(buffer_cinemri,3)
    tmp_imgs(:,:,i)=intshft(buffer_cinemri(:,:,i),mr_sh(i,:));
end

cinemri1=tmp_imgs(:,:,:);            %%%%%DEV TEST MARCH 18

 temp=strcat('Output/',int2str(slice),'/cinemri1.mat');
save(temp,'cinemri1');
