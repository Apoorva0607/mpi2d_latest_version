
function out=DR_core_new_affine(filename_ref,filename_tar,filename_out,filename_mask,sx,sy,step,a1,a2)

% increase Gauss wt and reduce SyN wt to avoid funny results

filename_out_nii=strcat(filename_out,'.nii');
%%%%%%%%%%%%%%%%%%%%%% OLD ANTS %%%%%%%%%%%%%%
 %eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/ANTS 2 -m CC\[',filename_ref,',',filename_tar,',1,3\] -i 100x100x10 -o ',filename_out_nii,' -t SyN\[0.1\] -r Gauss\[10,2\] -x ',filename_mask);
%eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/ANTS 2 -m CC\[',filename_ref,',',filename_tar,',1,3\] -i 100x100x10 -o ',filename_out_nii,' -t SyN\[4\] -r Gauss\[0,2\] -x ');
%eval(eval_command);
% 
% %  eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/Source/Binaries/bin/ANTS 2 -m MI\[',filename_ref,',',filename_tar,',1,15\] -i 0 -o ',filename_out_nii,' --rigid-affine true');
% %  eval(eval_command)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/092514/antsbin/bin/antsRegistration -d 2 -m MI[',filename_ref,',',filename_tar,',1,15] -o ',filename_out,' -t affine[0.1]');
 %eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/092514/antsbin/bin/antsRegistrationSyNQuick.sh -d 2 -f ',filename_ref,' -m ',filename_tar,' -o ',filename_out);
eval_command=horzcat('!/v/raid1a/gadluru/softs/ANTS/092514/antsbin/bin/antsRegistration -d 2 -r \[',filename_ref,' , ',filename_tar,', 1\]',...
                                                                                            ' -t rigid\[0.1\]',...
                                                                                            ' -m GC\[',filename_ref,',',filename_tar,', 1, 4]',...
                                                                                            ' -x ',filename_mask,...
                                                                                            ' -c \[100x50x10, 1e-30, 4\]',...
                                                                                            ' -s 3x2x1vox',...
                                                                                            ' -f 3x2x1',...
                                                                                            ' -o ',filename_out);

 eval(eval_command);

filename_out_warp=strcat(filename_out,'1Warp.nii.gz');
filename_out_affine=strcat(filename_out,'0GenericAffine.mat');
filename_tar_warp_mhd=strcat(filename_tar(1:end-4),'_warp.mhd');

eval_command2=horzcat('!/v/raid1a/gadluru/softs/ANTS/092514/antsbin/bin/antsApplyTransforms -d 2 -i ',filename_tar,' -o ',filename_tar_warp_mhd,' -r ',filename_ref,'  -n BSpline -t ',filename_out_affine);
eval(eval_command2);

  filename_tar_warp_raw=strcat(filename_tar(1:end-4),'_warp.raw');
 fid = fopen (filename_tar_warp_raw, 'r', 'ieee-le');
  temp_var=fread (fid,'double');
%  
 fclose (fid);
 out=reshape(temp_var,[sx sy]);
%  imagesc(out)
%  colormap gray
%out=temp_var;
return;


