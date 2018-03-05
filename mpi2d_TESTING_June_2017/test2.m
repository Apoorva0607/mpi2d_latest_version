function test2       

sliceNum=2;
temp=strcat('Output/',int2str(sliceNum),'/cinemri_curv_fit1.mat');
load(temp)
        for ifr=1:size(cinemri_curv_fit1,3)
            imagesc(cinemri_curv_fit1(:,:,ifr)),colormap gray,axis off
            pause(0.1)
        end