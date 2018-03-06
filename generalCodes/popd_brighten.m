function popd(images,beta, scaleFlag)

figure; colormap gray
 brighten(beta)
if exist('scaleFlag')
    maxx=max((images(:)));
    clim=[0 maxx]

for ii=1:size(images,3)
    imagesc(abs(images(:,:,ii)),clim); colorbar
    title(strcat('Scaled from 0 to max of all frames.  Frame/slice  ',int2str(ii)))
    pause
end
else


    while 5
for ii=1:size(images,3)
    imagesc(abs(images(:,:,ii)))
    title(int2str(ii))
    pause
end
    end
end
