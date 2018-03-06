
function out=DR_core(filename_ref,filename_tar,filename_out,sx,sy)

% increase Gauss wt and reduce SyN wt to avoid funny results

filename_out_nii=strcat(filename_out,'.nii');

 eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/ANTS 2 -m CC\[',filename_ref,',',filename_tar,',1,3\] -i 100x100x10 -o ',filename_out_nii,' -t SyN\[0.15\] -r Gauss\[1,2\]');
%eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/ANTS 2 -m CC\[',filename_ref,',',filename_tar,',1,3\] -i 100x100x10 -o ',filename_out_nii,' -t SyN\[4\] -r Gauss\[0,2\] -x ');
eval(eval_command);
% 
% %  eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/ANTS 2 -m MI\[',filename_ref,',',filename_tar,',1,15\] -i 0 -o ',filename_out_nii,' --rigid-affine true');
% %  eval(eval_command)
% '
filename_out_warp=strcat(filename_out,'Warp.nii');
filename_out_affine=strcat(filename_out,'Affine.txt');
filename_tar_warp_mhd=strcat(filename_tar(1:end-4),'_warp.mhd');


% eval_command2=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/WarpImageMultiTransform 2 ',filename_tar,' ',filename_tar_warp_mhd,' ',filename_out_affine,' -R ',filename_ref);
% eval(eval_command2)

eval_command2=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/WarpImageMultiTransform 2 ',filename_tar,' ',filename_tar_warp_mhd,' -R ',filename_ref,' --use-NN ',filename_out_warp,' ',filename_out_affine);
eval(eval_command2);


% eval_command2=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/WarpImageMultiTransform 2 ',filename_tar_warp_mhd,' ',filename_tar_warp_mhd,' -R ',filename_ref_fixed,' -i ',filename_out_warp_fixed,' ',filename_out_affine_fixed);
% eval(eval_command2);




filename_tar_warp_raw=strcat(filename_tar(1:end-4),'_warp.raw');
fid = fopen (filename_tar_warp_raw, 'r', 'ieee-le');
temp_var=fread (fid,'float');
fclose (fid);
out=reshape(temp_var,[sx sy]);
% imagesc(out)
% colormap gray
return;


