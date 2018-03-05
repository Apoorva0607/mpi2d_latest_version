clear all
clc
close all
warning off
oldoptions=optimset('fmincon');
options=optimset(oldoptions,'TolFun',0.00000001,'Display','off') ;
A=[];
B=[];
nonlcon=[];
x0(1)=double(700.0);
LB(1)=double(0);
UB(1)=double(2000.0);

list_dir=dir('CV*');

for dir_count=1:length(list_dir)

    cd(list_dir(dir_count).name);
    list_dicoms=dir('*.dcm');
    
    for dicom_count=1:length(list_dicoms)
        
        t1_series(:,:,dir_count,dicom_count)=double(dicomread(list_dicoms(dicom_count).name));
        
        
    end
    cd ..
end

SRT=[202,402];

for slice=1:size(t1_series,4)
slice
    t1_series_tmp=double(squeeze(t1_series(:,:,:,slice)));
    
    
for i=1:size(t1_series,1)
    
    parfor j=1:size(t1_series,2)
        A_scalar=t1_series_tmp(i,j,1);
        curve=squeeze(t1_series_tmp(i,j,2:3));
        
        if A_scalar>1
            t1_map(i,j,slice)= fmincon('T1_3param',x0,A,B,A,B,LB,UB,nonlcon,options,curve,A_scalar,SRT);
        else
            t1_map(i,j,slice)=0;
        end
        
    end
end

end

save('t1_map.mat','t1_map')
