function movie_play1(cinemri1)


figure

upsampleFactor=1;
for j=1:100
    for k=1:8 
    for i=1:size(cinemri1,4)
 i        
     
% a=interp2(cinemri1(:,:,i),upsampleFactor-1,'cubic');
a=cinemri1(:,:,k,i);
        imagesc((a(:,:)));axis image;
      colormap gray
axis off    
colorbar        
        
pause
    end   
    end
end
        