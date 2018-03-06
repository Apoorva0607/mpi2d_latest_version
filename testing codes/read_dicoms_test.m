% clear all
close all
clc
dicoms = dir('*.dcm');
h=fspecial('gaussian');
mask=zeros(144,144);
mask(70:80,70:80)=1;
 for j=1:length(dicoms)
    perf(:,:,j) = (double(dicomread(dicoms(j).name)));

 end
% tmp_pro=double((a1(:,:,1)+a1(:,:,2)+a1(:,:,3)+a1(:,:,4))./4);
% tmp_scale=mean(mean(tmp_pro));
% 
% for i=1:236
% a3(:,:,i)=(double(a1(:,:,i))./tmp_pro).*tmp_scale;
% end
% 
%  tmp=double((a3(:,:,9)+a3(:,:,10)+a3(:,:,8)+a3(:,:,11))./4);
%  
%  for j=1:length(dicoms)
%     a2(:,:,j) = double(a3(:,:,j))-tmp;
%  end
%  
%  parfor j=1:230
%  bw(:,:,j)=activecontour(a2(:,:,j),mask);
%  end