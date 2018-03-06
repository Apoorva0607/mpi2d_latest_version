% run in parallel with multiple cores

function [Registered_imgs]=ANTS_rigid(self_gated_imgs,model_imgs)



model_imgs=(model_imgs(:,:,1:end)/(max(max(max(model_imgs)))))*500;

scale_factor=500/max(max(max(self_gated_imgs)));
self_gated_imgs=self_gated_imgs*scale_factor;

sx=size(self_gated_imgs,1);
sy=size(self_gated_imgs,2);

rangex=1:sx;
rangey=1:sy;

nof=size(self_gated_imgs,3);

% [RV LV]=FindLVRV(self_gated_imgs,0);
% x=floor(mean([RV(1),LV(1)]));
% y=floor(mean([RV(2),LV(2)]));
% mask=zeros(sx,sy);
% mask(x-30:x+30,y-30:y+30)=1;
% 
 registered=zeros(size(self_gated_imgs));
% 
% filename_ref_mask=strcat('mask.raw');
%     fid = fopen (filename_ref_mask, 'wb', 'ieee-le');
%     fwrite (fid, mask, 'float');
%     fclose(fid);
%     
%   filename_mask=strcat('mask.mhd');
%     fid = fopen (filename_mask, 'wb');
%     header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_ref_mask);
%     fprintf(fid,'%s',header_info);
%     fclose(fid);  

parfor_progress(nof);
parfor i=1:nof

    ref_img=(model_imgs(:,:,i)/5);    % reference image is the fixed image
    %ref_img=(model_imgs./5);
    filename_ref_raw=strcat('ref_',int2str(i),'.raw');
    fid = fopen (filename_ref_raw, 'wb', 'ieee-le');
    fwrite (fid, ref_img(rangex,rangey), 'float');
    fclose(fid);
    
    filename_ref=strcat('ref_',int2str(i),'.mhd');
    fid = fopen (filename_ref, 'wb');
    header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_ref_raw);
    fprintf(fid,'%s',header_info);
    fclose(fid);
    
    tar_img=(self_gated_imgs(:,:,i)/5); % target image is the moving image
    filename_tar_raw=strcat('tar_',int2str(i),'.raw');
    fid = fopen (filename_tar_raw, 'wb', 'ieee-le');
    fwrite (fid, tar_img(rangex,rangey), 'float');
    fclose(fid);
    
    filename_tar=strcat('tar_',int2str(i),'.mhd');
    fid = fopen (filename_tar, 'wb');
    header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_tar_raw);
    fprintf(fid,'%s',header_info);
    fclose(fid);
        
     filename_out=strcat('out_',int2str(i));
%     filename_out_fixed=strcat('out_',int2str(ref_frame));
%     filename_ref_fixed=strcat('ref_',int2str(ref_frame),'.mhd');
    
    registered(:,:,i)=DR_core_rigid(filename_ref,filename_tar,filename_out,sx,sy);    
    parfor_progress;
end
parfor_progress(0);
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
% %     filename_out_fixed=strcat('out_',int2str(ref_frame));
% %     filename_ref_fixed=strcat('ref_',int2str(ref_frame),'.mhd');
%     
%     registered(:,:,i)=DR_core_invert_rigid(filename_ref,filename_tar,filename_out,sx,sy,filename_mask,filename_out_fixed,filename_ref_fixed);    
%    
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Registered_imgs=registered*5/scale_factor;

return;

