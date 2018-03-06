function [res]=mplot(xyt,time)
 
[sx,sy,sz]=size(xyt);
 
% xyt=permute(xyt,[1 2 3]);
 
 
 
colormap(gray)
 
 
tt=zeros(sx,sy);%
for i=1:sz,tt=xyt(:,:,i);imagesc((rot90(tt,-1)),[0 1.0]),axis image,
%     title(i),...%,[0 time]%,[0 3]
        set(gca,'XTick',[]),set(gca,'YTick',[]),...
        axis image,%title(i),
    
%    colorbar
    colormap(gray),
    %pause%(time),
    M(i)=getframe,
end
 
movie(M)
 
movie2avi(M,'movie0.avi','compression','Cinepak')
