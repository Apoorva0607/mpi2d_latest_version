close all;

figure(1); clf;

%load sens_comb_full_ro.mat

% modified EVRD 5/1/12
tmpFolders='/v/raid1/gadluru/of_for_by_people/Alexis/To_Alexis_for_paper/testing_data_new/CV_Radial7Off_flexICE_rcECG_pd_ILOT35_goldenratio_ungated_STRESS(*200)'

for ii=1:length(tmpFolders)
    
   dicomFiles = dir([tmpFolders(ii).name '/*.dcm']);
   for jj=1:length(dicomFiles)/2           
    temp(ii).images(:,:,jj) = dicomread(dicomFiles(jj).name]);
   end
end

    
    %process all folders in FoldersToDo
sens_comb_full=round(sens_comb_full*100);

z=imagesc(abs(sens_comb_full(:,:,1))),colormap gray,axis off,axis image,brighten(0.2);
% montage can also be used

set(gcf, 'Color' ,'w'); % set background to be white
f = getframe(figure(1)); % you get the whole picture including the axis
[im,map] = rgb2ind(f.cdata,1024);

% [im,map] = rgb2ind(f.cdata,1024,'nodither');
% [im,map] = gray2ind(f.cdata,512);

for k = 1:43
    k
    
    imagesc(abs(sens_comb_full(:,:,k))),colormap gray,axis off,axis image,brighten(0.2)
    set(gcf, 'Color' ,'w');
    axis tight
    set(gca,'nextplot','replacechildren','visible','off')
    f = getframe(figure(1)); % you get the whole picture including the axis
    im(:,:,1,k) = rgb2ind(f.cdata,map);
    
    %   im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
    %   im(:,:,1,k) = gray2ind(f.cdata,map);
    pause(0.1)
end

imwrite(im,map,'Output.gif','DelayTime',0.1)

