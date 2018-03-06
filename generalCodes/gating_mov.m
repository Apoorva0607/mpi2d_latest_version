clear all
close all
clc


vidobj=VideoWriter('/v/raid1a/dlikhite/selfGating.avi');
vidobj.FrameRate=10;
vidobj.Quality=100;
open(vidobj)
try
load('a.mat');
catch
dicoms = dir('*.dcm');
        if(isempty(dicoms))
            disp('No dicoms found.  Cannot do Temporal Smoothing registration');
            keyboard;
            return;
        end
        for t=1:size(dicoms,1)
            a(:,:,t) = mpi_upsampleImages(dicomread([dicoms(t).name]));
        end
%         save('a.mat','a');
end
a1=a(:,:,6:end);

% xc=128;
% yc=128;
% r1=40;
% r2=20;
% 
% 
% a=zeros(256,256,20);
% 
% for t=1:20
%  t   
% if mod(t,2)==0
%     
% for i=1:256
%     for j=1:256
%         
%         
%         if (i-xc).^2+(j-yc).^2<=r1.^2
%             a(i,j,t)=1;
%         end
%                
%         
%     end
% end
% 
% else
%  for i=1:256
%     for j=1:256
%         
%         
%         if (i-xc).^2+(j-yc).^2<=r2.^2
%             a(i,j,t)=1;
%         end
%                
%         
%     end
%  end
% end
%  
% end


for t=1:100
   h=subplot(2,1,1);
   imagesc(a1(:,:,t));
   hold on
   rectangle('Position',[30,50,80,60],'EdgeColor','r')
   axis square
   
 set(gca,'ytick',[]) ;
    set(gca,'xtick',[]) ; 
   colormap gray
   brighten(0.3)
   subplot(2,1,2)
    curv(t)=sum(sum(a1(50:110,30:110,t)));
    plot(curv,'-o','MarkerEdgeColor','r','MarkerFaceColor','r')
%    axis square
   xlim([1 100]);
   ylim([3.5e5 8e5])
    
%     ylim([1000 5500])
    % xlim([10 100])
    
    frame1=getframe(gcf);
    writeVideo(vidobj,frame1);
   
end
 close(vidobj)           