function movie_play(cinemri1)


figure

upsampleFactor=1;
for j=1:100
    for i=1:size(cinemri1,3)
 i        
      
% a=interp2(cinemri1(:,:,i),upsampleFactor-1,'cubic');
a=cinemri1(:,:,i);
        imagesc((a(:,:)));axis image;
      colormap gray
axis off    
colorbar        
        
pause
        
    end
end
        