
clear

% load /v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P102915/ReconData/mat_files/hybrid/MID_92/dia/more_rays/vbm4d/new/registered_imgs.mat
% load /v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P102915/ReconData/mat_files/hybrid/MID_92/dia/more_rays/vbm4d/new/add_tTV/set3/registered_imgs.mat
% dia1=Registered_imgs1;
% dia2=Registered_imgs2;
% dia3=Registered_imgs3;

% load /v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P102915/ReconData/mat_files/hybrid/MID_92/sys/more_rays/vbm4d/new/registered_imgs.mat
% load /v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P102915/ReconData/mat_files/hybrid/MID_92/sys/more_rays/vbm4d/new/add_tTv/set3/registered_imgs.mat
% sys1=Registered_imgs1;
% sys2=Registered_imgs2;
% sys3=Registered_imgs3;

sx=288;
load /v/raid1a/dlikhite/collabs/Test/Grant' registration for movie'/op_img_1000_1_Droid_new.mat
dia1=op_img(sx/4+1:3*sx/4,sx/4+1:3*sx/4,:);
load /v/raid1a/dlikhite/collabs/Test/Grant' registration for movie'/op_img_1000_2_Droid_new.mat
dia2=op_img(sx/4+1:3*sx/4,sx/4+1:3*sx/4,4:end);
load /v/raid1a/dlikhite/collabs/Test/Grant' registration for movie'/op_img_1000_3_Droid_new.mat
dia3=op_img(sx/4+1:3*sx/4,sx/4+1:3*sx/4,4:end);

load /v/raid1a/dlikhite/collabs/Test/Grant' registration for movie'/op_img_4000_1_Droid_new.mat
sys1=op_img(sx/4+1:3*sx/4,sx/4+1:3*sx/4,4:end);
load /v/raid1a/dlikhite/collabs/Test/Grant' registration for movie'/op_img_4000_2_Droid_new.mat
sys2=op_img(sx/4+1:3*sx/4,sx/4+1:3*sx/4,4:end);
load /v/raid1a/dlikhite/collabs/Test/Grant' registration for movie'/op_img_4000_3_Droid_new.mat
sys3=op_img(sx/4+1:3*sx/4,sx/4+1:3*sx/4,4:end);

clear temp M
iptsetpref('ImshowBorder','tight'); % use this to get rid of white border for each frame

for i=1:60
    temp(:,:,1,1)=rot90(sys1(16:end-15,16:end-15,i));
    temp(:,:,1,2)=rot90(dia1(16:end-15,16:end-15,i));
    temp(:,:,1,3)=rot90(sys2(16:end-15,16:end-15,i));
    temp(:,:,1,4)=rot90(dia2(16:end-15,16:end-15,i));
    temp(:,:,1,5)=rot90(sys3(16:end-15,16:end-15,i));
    temp(:,:,1,6)=rot90(dia3(16:end-15,16:end-15,i));
    temp=temp/max(temp(:));
%     set(gca,'position',[0 0 1 1],'units','normalized')
    
    montage(temp,'size',[3 2]),brighten(0.4)
%     pause
    M(i)=getframe(gcf);
end
movie2avi(M,'GM1_newRegi_Droid.avi','fps',10)

% keyboard
% % load /v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P102915/ReconData/mat_files/hybrid/MID_92/dia/more_rays/vbm4d/new/new_reduced_non_k_space_sms_MID_92_set_1iter_no100.mat
% load /v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P102915/ReconData/mat_files/hybrid/MID_92/dia/more_rays/vbm4d/new/add_tTV/set2/new_reduced_non_k_space_sms_MID_92_set_1.mat
% sx=288;
% dia1=abs(new_reduced_non_k_space_slice_1(sx/4+1:3*sx/4,sx/4+1:3*sx/4,:));
% dia2=abs(new_reduced_non_k_space_slice_2(sx/4+1:3*sx/4,sx/4+1:3*sx/4,:));
% dia3=abs(new_reduced_non_k_space_slice_3(sx/4+1:3*sx/4,sx/4+1:3*sx/4,:));
% 
% % load /v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P102915/ReconData/mat_files/hybrid/MID_92/sys/more_rays/vbm4d/new/new_reduced_non_k_space_sms_MID_92_set_1iter_no100.mat
% load /v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P102915/ReconData/mat_files/hybrid/MID_92/sys/more_rays/vbm4d/new/add_tTv/set2/new_reduced_non_k_space_sms_MID_92_set_1.mat
% sys1=abs(new_reduced_non_k_space_slice_1(sx/4+1:3*sx/4,sx/4+1:3*sx/4,:));
% sys2=abs(new_reduced_non_k_space_slice_2(sx/4+1:3*sx/4,sx/4+1:3*sx/4,:));
% sys3=abs(new_reduced_non_k_space_slice_3(sx/4+1:3*sx/4,sx/4+1:3*sx/4,:));
% 
% for i=1:60
%     temp(:,:,1,1)=rot90(sys1(16:end-25,26:end-20,i));
%     temp(:,:,1,2)=rot90(dia1(16:end-25,26:end-20,i));
%     temp(:,:,1,3)=rot90(sys2(16:end-25,26:end-20,i));
%     temp(:,:,1,4)=rot90(dia2(16:end-25,26:end-20,i));
%     temp(:,:,1,5)=rot90(sys3(16:end-25,26:end-20,i));
%     temp(:,:,1,6)=rot90(dia3(16:end-25,26:end-20,i));
%     temp=temp/max(temp(:));
%     montage(temp,'size',[3 2]),brighten(0.4)    
%     M(i)=getframe(gcf);
% end
% movie2avi(M,'P102916_new3.avi','fps',10)



% % new video writer
% 
% vidObj = VideoWriter('peaks.avi','Uncompressed AVI');
% vidObj.FrameRate=10;
% open(vidObj);
% for i=1:60
%     temp(:,:,1,1)=rot90(sys1(16:end-15,16:end-15,i));
%     temp(:,:,1,2)=rot90(dia1(16:end-15,16:end-15,i));
%     temp(:,:,1,3)=rot90(sys2(16:end-15,16:end-15,i));
%     temp(:,:,1,4)=rot90(dia2(16:end-15,16:end-15,i));
%     temp(:,:,1,5)=rot90(sys3(16:end-15,16:end-15,i));
%     temp(:,:,1,6)=rot90(dia3(16:end-15,16:end-15,i));
%     temp=temp/max(temp(:));
%     montage(temp,'size',[3 2]),brighten(0.4)
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
% end
% 
