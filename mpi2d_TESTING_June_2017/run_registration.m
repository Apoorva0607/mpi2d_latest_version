
function run_registration(fullcinemri1,rangey,rangex,referenceFrame, numSkip,numPreContrast_bld,numPreContrast_tiss,sliceNum)
num_slices=1;

% model based registration

% Input: buffer_cinemri.mat (full FOV image) in the Output folder
% Output: Output/cinemri1.mat (registered movie)

% set registration parameters here
numSkip=0;

noi=1;  % total number of MBR iterations


    
        
    try
        temp=strcat('Output/',int2str(sliceNum),'/cinemri_curv_fit1.mat');
        load(temp)
%         for ifr=1:size(cinemri1,3)
%             imagesc(cinemri1(:,:,ifr)),colormap gray,axis off
%             pause(0.1)
%         end
     catch
        upsample_flag=0;upsample_factor=2;
%         figure(1);imagesc(fullcinemri1(:,:,20));colormap gray;
%         mask=roipoly;
%         nof=size(fullcinemri1,3);
%         for i=1:nof
%             fullcinemri(:,:,i)=mask.*fullcinemri1(:,:,i);
%         end
        MBR_main(fullcinemri1,rangex,rangey,referenceFrame,numSkip,numPreContrast_bld,numPreContrast_tiss,noi,sliceNum,upsample_flag,upsample_factor);
       end    
%         


