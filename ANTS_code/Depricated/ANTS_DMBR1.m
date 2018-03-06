% run in parallel with multiple cores

function [Registered_imgs]=ANTS_DMBR1(self_gated_imgs,model_imgs,mask,step,a1,a2)

%ref_frame=10;

% model_imgs=(model_imgs(:,:,1:end)/(max(max(max(model_imgs)))))*500;
% 
% scale_factor=500/max(max(max(self_gated_imgs)));
% self_gated_imgs=self_gated_imgs*scale_factor;

sx=size(self_gated_imgs,1);
sy=size(self_gated_imgs,2);

rangex=1:sx;
rangey=1:sy;

nof=size(self_gated_imgs,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [RV LV]=FindLVRV(self_gated_imgs,0);
% x=floor(mean([RV(1),LV(1)]));
% y=floor(mean([RV(2),LV(2)]));
% mask=zeros(sx,sy);
% mask(x-20:x+20,y-20:y+20)=1;
%imagesc(self_gated_imgs(:,:,ref_frame));colormap gray;axis image
% mask=roipoly;

% [RV,LV] = FindLVRV(self_gated_imgs,0);  
% 
%   % just center on RV: 
%   
% %    sy=size(self_gated_imgs,1);
% %    sx=size(self_gated_imgs,2);
% %   RV=[round(sx/2),round(sy/2)];
% %  LV=[round(sx/2),round(sy/2)];
%              windowWidth = sx/1.8;
%              windowHeight = sy/2;
%              H_cen(2)=(RV(2)+LV(2))/2;
%              H_cen(1)=(RV(1)+LV(1))/2;
%              rangex1 = round(H_cen(2) - windowWidth/2):round(H_cen(2) + windowWidth/2);
%              rangey1 = round(H_cen(1) - windowHeight/2):round(H_cen(1) + windowHeight/2);
% % rangex=30:sx-30;
% % rangey=30:sy-30;
% % 
% 
% rangex2=rangex1;
% rangey2=rangey1;
%             sz_adjust=1
% %             rangex2 = round(RV(1) - sz_adjust*windowWidth):round(RV(1) + sz_adjust*windowWidth/2);
%             if round(H_cen(2) - sz_adjust*windowWidth/2) < 5
%                rangex2=5:round(H_cen(2) + sz_adjust*windowHeight/2);
%             end
%             if round(H_cen(2) + sz_adjust*windowWidth/2) >sx
%                rangex2=round(H_cen(2) - sz_adjust*windowWidth/2):sx-5;
%             end
%             
% %            rangey2 = round(RV(2) - sz_adjust*windowHeight/2):round(RV(2) + sz_adjust*windowHeight);
%             if round(H_cen(1) - sz_adjust*windowHeight/2) < 5
%                rangey2=5:round(H_cen(1) + sz_adjust*windowHeight/2);
%             end
%             if round(H_cen(1) + sz_adjust*windowHeight/2)>sy
%                rangey2=round(H_cen(1) - sz_adjust*windowHeight/2):sy-5;
%             end
% 
% rangex1=rangex2;
% rangey1=rangey2;
% mask=zeros(sx,sy);
% mask(rangey1,rangex1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
registered=zeros(size(self_gated_imgs));

% filename_ref_mask=strcat('/working/ANTS_tmp/mask.raw');
%     fid = fopen (filename_ref_mask, 'wb', 'ieee-le');
%     fwrite (fid, mask, 'float');
%     fclose(fid);
%     
%   filename_mask=strcat('/working/ANTS_tmp/mask.mhd');
%     fid = fopen (filename_mask, 'wb');
%     header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_ref_mask);
%     fprintf(fid,'%s',header_info);
%     fclose(fid);  

%parfor_progress(nof);

%registered(:,:,1)=self_gated_imgs(:,:,1)/5;
parfor i=1:nof

    ref_img=(model_imgs(:,:,i)/5);    % reference image is the fixed image
    %ref_img=(registered(:,:,i-1));
    filename_ref_raw=strcat('/v/raid1a/apedgaon/ANTS/ANTS_tmp/ref_',int2str(i),'.raw');
    fid = fopen (filename_ref_raw, 'wb', 'ieee-le');
    fwrite (fid, ref_img(rangex,rangey), 'float');
    fclose(fid);
    
    filename_ref=strcat('/v/raid1a/apedgaon/ANTS/ANTS_tmp/ref_',int2str(i),'.mhd');
    fid = fopen (filename_ref, 'wb');
    header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_ref_raw);
    fprintf(fid,'%s',header_info);
    fclose(fid);
    
    tar_img=(self_gated_imgs(:,:,i)/5); % target image is the moving image
    filename_tar_raw=strcat('/v/raid1a/apedgaon/ANTS/ANTS_tmp/tar_',int2str(i),'.raw');
    fid = fopen (filename_tar_raw, 'wb', 'ieee-le');
    fwrite (fid, tar_img(rangex,rangey), 'float');
    fclose(fid);
    
    filename_tar=strcat('/v/raid1a/apedgaon/ANTS/ANTS_tmp/tar_',int2str(i),'.mhd');
    fid = fopen (filename_tar, 'wb');
    header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_tar_raw);
    fprintf(fid,'%s',header_info);
    fclose(fid);
        
    filename_out=strcat('/v/raid1a/apedgaon/ANTS/ANTS_tmp/out_',int2str(i));
%     filename_out_fixed=strcat('out_',int2str(ref_frame));
%     filename_ref_fixed=strcat('ref_',int2str(ref_frame),'.mhd');
   
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%% MASK %%%%%%%%%%%%%%%%%
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     mask_img=(mask); % target image is the moving image
    filename_mask_raw=strcat('/v/raid1a/apedgaon/ANTS/ANTS_tmp/mask_',int2str(i),'.raw');
    fid = fopen (filename_mask_raw, 'wb', 'ieee-le');
    fwrite (fid, mask_img(rangex,rangey), 'float');
    fclose(fid);
    
    filename_mask=strcat('/v/raid1a/apedgaon/ANTS/ANTS_tmp/mask_',int2str(i),'.mhd');
    fid = fopen (filename_mask, 'wb');
    header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_mask_raw);
    fprintf(fid,'%s',header_info);
    fclose(fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    registered(:,:,i)=DR_core_new1(filename_ref,filename_tar,filename_out,filename_mask,sx,sy,step,a1,a2);    
    %parfor_progress;
end
%parfor_progress(0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parfor i=1:nof
% 
%     ref_img=(model_imgs(:,:,i)/5);    % reference image is the fixed image
%     %ref_img=(model_imgs./5);
%     filename_ref_raw=strcat('ref_',int2str(i),'.raw');
%     fid = fopen (filename_ref_raw, 'wb', 'ieee-le');
%     fwrite (fid, ref_img(rangex,rangey), 'float');
%     fclose(fid);
%     
%     filename_ref=strcat('ref_',int2str(i),'.mhd');
%     fid = fopen (filename_ref, 'wb');
%     header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_ref_raw);
%     fprintf(fid,'%s',header_info);
%     fclose(fid);
%     
%     tar_img=(self_gated_imgs(:,:,i)/5); % target image is the moving image
%     filename_tar_raw=strcat('tar_',int2str(i),'.raw');
%     fid = fopen (filename_tar_raw, 'wb', 'ieee-le');
%     fwrite (fid, tar_img(rangex,rangey), 'float');
%     fclose(fid);
%     
%     filename_tar=strcat('tar_',int2str(i),'.mhd');
%     fid = fopen (filename_tar, 'wb');
%     header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_tar_raw);
%     fprintf(fid,'%s',header_info);
%     fclose(fid);
%         
%     filename_out=strcat('out_',int2str(i));
%     filename_out_fixed=strcat('out_',int2str(ref_frame));
%     filename_ref_fixed=strcat('ref_',int2str(ref_frame),'.mhd');
%     
%     registered(:,:,i)=DR_core_invert(filename_ref,filename_tar,filename_out,sx,sy,filename_mask,filename_out_fixed,filename_ref_fixed);    
%    
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Registered_imgs=registered*5;
%Registered_imgs=registered;
return;
