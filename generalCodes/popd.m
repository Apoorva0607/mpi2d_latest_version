function popd(images,titlet,scaleFlag )

colormap gray

if ~exist('titlet')
    titlet=' ';
end

if exist('scaleFlag')
    maxx=max((images(:)));
    clim=[0 maxx]
    clim=[0 maxx*0.7]

for ii=1:size(images,3)
    imagesc(abs(images(:,:,ii)),clim); colorbar
    title(strcat('Scaled from 0 to max of all frames.  Frame/slice  ',int2str(ii)))
    pause
end
else

for ii=1:size(images,3)
    imagesc(abs(images(:,:,ii)))
    title(strcat(titlet, int2str(ii)))
    pause
end
end
