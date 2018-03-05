clc;
close all;
clear all;
list=dir(pwd);
    folder=pwd;
    files = dir([folder '/*.par']);
    outpath='Output/';

for i=1:length(files)-2
    seriesNum=1000*i+1000;
    sliceNum=i+1;
    seriesNumAIF=100;
    sliceNumAIF=1;
    scaleAIF=1;
    filename=strcat(outpath,'est_curves.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
    load(filename);
     filename1=strcat(outpath,'tisscurve.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
    load(filename1);
    filename2=strcat(outpath,'ttimes',int2str(seriesNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
    load(filename2);
    [nRegs, nTimes]=size(tisscurve');
%     figure(i);
for k=1:nRegs
    h=figure(i+1); hold on; set(gcf,'position',[100 100 800 400])
    titlet='slice ';
    title(strcat(titlet, int2str(i+1)))
   
    ColOrd = get(gca,'ColorOrder');
    plot(ttimes,est_curves(1:nTimes,k),'Color',ColOrd(k,:),'LineWidth',2)
    plot(ttimes,tisscurve(:,k),'x--','Color',ColOrd(k,:),'LineWidth',1.3)
    ylim([0 1.5]);
end
end