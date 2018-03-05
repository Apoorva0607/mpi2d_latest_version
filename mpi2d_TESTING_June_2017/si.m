clc
close all
clear all
infile='/v/raid1a/group/MRIdata/Cardiac/Prisma/P082516/RegisteredData/RegisteredData_P082516_MID57.mat';
% infile='/v/raid1a/ytian/MRIdata/Cardiac/3D_Data_For_Jason/P082516/processing/SET08812_MID057_GROG_Perfusion_AIF_3D_echo_1_170818_160051.mat';
load(infile);
% for i=1:size(Image,3)
% I(:,:,i)=rot90(Image(:,:,i),2);
% end

%infile='/v/raid1a/ytian/MRIdata/Cardiac/3D_Data_For_Jason/P082516/processing/SET08813_MID057_GROG_Perfusion_AIF_3D_echo_2_170818_160234.mat';
% load(infile);
 Image=Data.RegisteredNormalAIF;
Image1=Image(:,:,:,1);
Image2=Image(:,:,:,2);
I=Image1;
J=Image2;
% for j=1:size(Image,3)
% J(:,:,j)=rot90(Image(:,:,j),2);
% end

% infile='/v/raid1a/ytian/MRIdata/Cardiac/3D_Data_For_Jason/P082516/processing/SET08813_MID057_GROG_Perfusion_AIF_3D_echo_2_170818_160234.mat';
% load(infile);
% J=Image2;
I1=I(:,:,32);
J1=J(:,:,32);
outpath='/v/raid1a/apedgaon/MRIdata/Cardiac/Prisma_version1/P082516_reg/Processing_rest/Output/';
%%% if the ROIs need to be drawn separately on the image the flag man
%%% change to 1
man=0;
if(man==0)
seriesNum=100;
sliceNum=1;
LV = load([outpath 'blood_polyCoords.study' int2str(seriesNum) '.slice' int2str(sliceNum)]); %%%% already drawn contours for AIF and saved by stage 3.11
bw = mpi_roipoly(I1,LV(:,1),LV(:,2));
figure(1);imagesc(I1);colormap gray;
hold on
plot(LV(:,1),LV(:,2));
title 'echo1 oldrecon'
figure(2);imagesc(J1);colormap gray;
hold on
plot(LV(:,1),LV(:,2));
title 'echo2 oldrecon'
else
  figure(1);imagesc(I1);colormap gray;
 bw=roipoly;
end

for i=1:size(I,3)
    temp1=I(:,:,i);
    temp2=J(:,:,i);
    echo1(i)=mean(temp1(bw>0));
    echo2(i)=mean(temp2(bw>0));
   
end
figure(3);plot(echo1,'R');
hold on
figure(3);plot(echo2,'B');
xlabel('time')
ylabel('SI')
if(man==0)
title 'comparison of two echos old recon'
else
title 'comparison of two echos without using mpi2d stg=3.11'
end