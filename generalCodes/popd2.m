function popd(images1,images2,scaleFlag)

colormap gray

if exist('scaleFlag')
    maxx=max((images1(:)));
    clim=[0 maxx]

for ii=1:size(images1,3)
    subplot(2,2,1)
    imagesc(abs(images1(:,:,ii)),clim); colorbar
    subplot(2,2,2)
    imagesc(abs(images2(:,:,ii)),clim); colorbar
    title(strcat('Scaled from 0 to max of all frames.  Frame/slice  ',int2str(ii)))
    pause
end
else

for ii=1:size(images1,3)
    %imagesc(abs(images(:,:,ii)))
    subplot(2,2,1)
    imagesc(abs(images1(:,:,ii))); colorbar
    subplot(2,2,2)
    imagesc(abs(images2(:,:,ii))); colorbar
    
    title(int2str(ii))
    pause
end
end
