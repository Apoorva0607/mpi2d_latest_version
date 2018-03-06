
function out=DR_core_invert(filename_ref,filename_tar,filename_out,sx,sy,filename_mask,filename_out_fixed,filename_ref_fixed)

% increase Gauss wt and reduce SyN wt to avoid funny results

filename_out_nii=strcat(filename_out,'.nii');

% eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/ANTS 2 -m CC\[',filename_ref,',',filename_tar,',1,2\] -i 100x100x10 -o ',filename_out_nii,' -t SyN\[0.1\] -r Gauss\[1,0\] -x ',filename_mask);
% eval(eval_command);
% 
% %  eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/ANTS 2 -m MI\[',filename_ref,',',filename_tar,',1,15\] -i 0 -o ',filename_out_nii,' --rigid-affine true');
% %  eval(eval_command)
% 
filename_out_warp=strcat(filename_out,'Warp.nii');
filename_out_affine=strcat(filename_out,'Affine.txt');
filename_tar_warp_mhd=strcat(filename_tar(1:end-4),'_warp.mhd');

filename_out_warp_fixed=strcat(filename_out_fixed,'Warp.nii');
filename_out_affine_fixed=strcat(filename_out_fixed,'Affine.txt');




% eval_command2=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/WarpImageMultiTransform 2 ',filename_tar,' ',filename_tar_warp_mhd,' ',filename_out_affine,' -R ',filename_ref);
% eval(eval_command2)

% eval_command2=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/WarpImageMultiTransform 2 ',filename_tar,' ',filename_tar_warp_mhd,' -R ',filename_ref,' --use-NN ',filename_out_warp,' ',filename_out_affine);
% eval(eval_command2);


eval_command2=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/WarpImageMultiTransform 2 ',filename_tar_warp_mhd,' ',filename_tar_warp_mhd,' -R ',filename_ref_fixed,' --use-NN -i ',filename_out_affine_fixed);
eval(eval_command2);

eval_command2=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/WarpImageMultiTransform 2 ',filename_tar_warp_mhd,' ',filename_tar_warp_mhd,' -R ',filename_ref_fixed,' --use-NN -i ',filename_out_warp_fixed);
eval(eval_command2);


filename_tar_warp_raw=strcat(filename_tar(1:end-4),'_warp.raw');
fid = fopen (filename_tar_warp_raw, 'r', 'ieee-le');
temp_var=fread (fid,'float');
fclose (fid);
out=reshape(temp_var,[sx sy]);
% imagesc(out)
% colormap gray
return;

