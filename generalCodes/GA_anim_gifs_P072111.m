clear 
figure(1); clf;

tmppwd=pwd

tmpFolders='/v/raid1/ed/MRIdata/Cardiac/Verio/P072111/ReconData/CV_Radial7Off_flexICE_rcECG_pdTEST_reconILOT35(2*)'
chdir('/v/raid1/ed/MRIdata/Cardiac/Verio/P072111/ReconData/')
ss=dir(tmpFolders)
ss=ss(1:5);
ss(1).name='CV_Radial7Off_flexICE_rcECG_pdTEST_reconILOT35(19000)'
ss(2).name='CV_Radial7Off_flexICE_rcECG_pdTEST_reconILOT35(20000)'
ss(3).name='CV_Radial7Off_flexICE_rcECG_pdTEST_reconILOT35(21000)'
ss(4).name='CV_Radial7Off_flexICE_rcECG_pdTEST_reconILOT35(22000)'
ss(5).name='CV_Radial7Off_flexICE_rcECG_pdTEST_reconILOT35(23000)'

startframe=50
for ii=1:length(ss)
    
   %dicomFiles = dir([ss(ii).name '/*.dcm']);
   tmpFiles = dir([ss(ii).name '/*.dcm']);
   for iframe=startframe:108  %length(tmpFiles)/2                 
    stress_images(:,:,iframe-startframe+1,ii) = dicomread([ss(ii).name '/' tmpFiles(iframe).name]);
   end
end


tmpFolders='/v/raid1/ed/MRIdata/Cardiac/Verio/P072111/ReconData/CV_Radial7Off_flexICE_rcECG_pdTEST_reconILOT35_8.5ml_ungated(5*)'
ss=dir(tmpFolders)
startframe=190
for ii=1:length(ss)
    ii
   tmpFiles = dir([ss(ii).name '/*.dcm']);
   for iframe=startframe:380 %1:length(tmpFiles)/2          
    rest_images(:,:,iframe-startframe+1,ii) = dicomread([ss(ii).name '/' tmpFiles(iframe).name]);
   end
end

chdir(tmppwd)

cropx=21:108
cropy=21:108
stress_images=5*stress_images(cropx,cropy,:,:);
rest_images=5*rest_images(cropx,cropy,:,:);

% all_dia=abs(all_dia)/8;
% all_dia=all_dia(21:108,21:108,:,:);


% all_sys=abs(all_sys)/8;
% all_sys=all_sys(21:108,21:108,:,:);
% 
% b(:,:,1,1)=sl(:,:,10);
% b(:,:,1,2)=sl(:,:,10);
% b(:,:,1,3)=sl(:,:,10);
% b(:,:,1,4)=sl(:,:,10);
% b(:,:,1,5)=sl(:,:,10);
% 
% b(:,:,1,6)=all_sys(:,:,10);
% b(:,:,1,7)=all_sys(:,:,10);
% b(:,:,1,8)=all_sys(:,:,10);
% b(:,:,1,9)=all_sys(:,:,10);
% b(:,:,1,10)=all_sys(:,:,10);
% 
% % b(:,:,1,11)=all_sys(:,:,10);
% % b(:,:,1,12)=all_sys(:,:,10);
% % b(:,:,1,13)=all_sys(:,:,10);
% % b(:,:,1,14)=all_sys(:,:,10);
% b(:,:,1,15)=all_sys(:,:,10);

    k=1 
    offset=18
    offset=0
    
%     
%     b(:,:,1,1)=stress_images(:,:,k+offset,3);
%     b(:,:,1,2)=stress_images(:,:,k+offset,5)/1.35;
%     b(:,:,1,3)=stress_images(:,:,k+offset,2);
%     b(:,:,1,4)=stress_images(:,:,k+offset,4)/1.35;
%     b(:,:,1,5)=stress_images(:,:,k+offset,1)*1.5;
%     
%     
%     b(:,:,1,6)=rest_images(:,:,k+offset,3);    
%     b(:,:,1,7)=rest_images(:,:,k+offset,5)/1.35;    
%     b(:,:,1,8)=rest_images(:,:,k+offset,2);    
%     b(:,:,1,9)=rest_images(:,:,k+offset,4)/1.35;    
%     b(:,:,1,10)=rest_images(:,:,k+offset,1);
%     
%     
% z=montage(b,'size',[2 5]);brighten(0.2)
% % montage can also be used
% 
% set(gcf, 'Color' ,'w'); % set background to be white
% f = getframe(figure(1)); % you get the whole picture including the axis
% [im,map] = rgb2ind(f.cdata,1024);
% 
% % [im,map] = rgb2ind(f.cdata,1024,'nodither');
% % [im,map] = gray2ind(f.cdata,512);
% 
% offset=18 %10+10+8;
% 
% map
% keyboard

for k = 1:159  %90
    k
        
%     b(:,:,1,1)=stress_images(:,:,k+offset,1)*1.3;
%     b(:,:,1,2)=stress_images(:,:,k+offset,2);
%     b(:,:,1,3)=stress_images(:,:,k+offset,3);
%     b(:,:,1,4)=stress_images(:,:,k+offset,4)/1.3;
%     b(:,:,1,5)=stress_images(:,:,k+offset,5)/1.44;
    
    b(:,:,1,1)=rest_images(:,:,k+offset,1)*1.1;
    b(:,:,1,2)=rest_images(:,:,k+offset,4)/1.3;
    b(:,:,1,3)=rest_images(:,:,k+offset,2);
    b(:,:,1,4)=rest_images(:,:,k+offset,5)/1.44;
    b(:,:,1,5)=rest_images(:,:,k+offset,3);
   
    
    
    %montage(b/9.8,'size',[2 5]);brighten(0.2)
    montage(1.6*b,'size',[1 5], 'displayRange', [0 5000]);brighten(0.3)
    mm=0:.01:1;
    mm=[mm;mm;mm]';
    set(gcf, 'Color' ,'w');
    axis tight
    set(gca,'nextplot','replacechildren','visible','off')
    f = getframe(figure(1)); % you get the whole picture including the axis
    im(:,:,1,k) = rgb2ind(f.cdata,mm);
    
    %   im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
    %   im(:,:,1,k) = gray2ind(f.cdata,map);
    pause(0.1)
end

%imwrite(im,mm,'ungated.gif','DelayTime',0.1)

