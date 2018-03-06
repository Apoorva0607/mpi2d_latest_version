function popdiff(img1,img2,scaleFlag)

colormap gray

% check if complex
images=img1-img2;
images=abs(images);

if exist('scaleFlag')
    maxx=max((images(:)));
    minn=min((images(:)));
    clim=[minn maxx]

for ii=1:size(images,3)
    imagesc((images(:,:,ii)),clim); colorbar
    title(strcat('Scaled from min to max of all frames.  Frame/slice  ',int2str(ii)))
    pause
end
else

for ii=1:size(images,3)
    imagesc((images(:,:,ii)))
    title(int2str(ii))
    pause
end
end
