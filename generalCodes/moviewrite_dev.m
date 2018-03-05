function moviewrite_dev(in_file)

vidobj=VideoWriter('regmovie2.mov');
vidobj.FrameRate=5;
vidobj.Quality=100;
open(vidobj)

in_file_up=mpi_upsampleImages(in_file,2);
for i=1:1:size(in_file_up,3)
   
    
    
%     imagesc((in_file_up(:,:,i)));
    imshow((in_file_up(:,:,i)),[0.02 0.5]);
    
      colormap gray
%     colorbar
%     axis image
%     axis off
    title(['Slice ',num2str(i)]);
%     brighten(0.4)
    frame1=getframe(gcf);
    size(frame1)
    writeVideo(vidobj,frame1);
end

close(vidobj)