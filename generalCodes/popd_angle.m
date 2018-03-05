function popd_angle(images,scaleFlag)

colormap jet   %gray

if exist('scaleFlag')
    maxx=max(angle(images(:)));
    minn=min(angle(images(:)));
    clim=[minn maxx]

for ii=1:size(images,3)
    imagesc(angle(images(:,:,ii)),clim); colorbar
    title(strcat('Scaled from min to max of all frames.  Frame/slice  ',int2str(ii)))
    pause
end
else

for ii=1:size(images,3)
    imagesc(angle(images(:,:,ii)))
    title(strcat('Phase ', int2str(ii)))
    pause
end
end
