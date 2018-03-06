% Model based registration

% Input: buffer_cinemri.mat (full FOV image) in the current folder
% Output: Output/cinemri1.mat (registered movie)


try
    matlabpool close
end
try
    matlabpool local 8
end

clear;

num_slices=3;

% model based registration

addpath MBR_Code




% set registration parameters here
numSkip=0;

noi=1;  % total number of MBR iterations

for i=1:num_slices
    
    filename=strcat('buffer_cinemri_slice',int2str(i),'.mat');
    load(filename)
    
    try
        flname=strcat('params_slice_',int2str(i),'.mat');
        load(flname)
        
    catch
        
        messg=strcat('Specify reference frame for slice#',int2str(i),' :');
        referenceFrame=input(messg);
        
        imagesc(buffer_cinemri(:,:,referenceFrame)),colormap gray,brighten(0.2),title(strcat('Slice # ',int2str(i))),axis image
        
        disp('Draw square region for cropping')
        
        [bw x y]=roipoly;
        rangex=round(min(x)):round(max(x));
        rangey=round(min(y)):round(max(y));
        
        disp('View Cine')
        
        for fr_no=1:size(buffer_cinemri,3)
            imagesc(buffer_cinemri(:,:,fr_no)),colormap gray,brighten(0.2),title(int2str(fr_no)),axis image
            pause
        end
        
        numPreContrast_bld=input('Specify number of pre-contrast frames for deltaSI curves :');
        numPreContrast_tiss=numPreContrast_bld;
        
        params.referenceFrame=referenceFrame;
        params.rangex=rangex;
        params.rangey=rangey;
        params.numPreContrast_bld=numPreContrast_bld;
        params.numPreContrast_tiss=numPreContrast_tiss;
        
        flname=strcat('params_slice_',int2str(i),'.mat');
        save(flname,'params');
    end
    
    try
        temp=strcat('Output/',int2str(i),'/cinemri1.mat');  % if registration already done
        load(temp)
        for ifr=1:size(cinemri1,3)
            imagesc(cinemri1(:,:,ifr)),colormap gray,axis off
            pause
        end
    catch
        MBR_main(params.rangex,params.rangey,params.referenceFrame,numSkip,params.numPreContrast_bld,params.numPreContrast_tiss,noi,i);
    end    
        
end

try
    matlabpool close
end