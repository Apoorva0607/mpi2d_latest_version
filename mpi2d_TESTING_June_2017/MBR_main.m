function MBR_main(buffer_cinemri,rangey,rangex,referenceFrame,numSkip,numPreContrast_bld,numPreContrast_tiss,noi,slice,upsample_flag,upsample_factor)

temp=strcat('Output/',int2str(slice));
mkdir(temp)


%%% stage 2.1
if(upsample_flag==1)
    
    if(upsample_factor==2)
        for i=1:size(buffer_cinemri,3)
            temp_buffer_cinemri(:,:,i)=interp2(buffer_cinemri(:,:,i),1,'cubic');
        end
        buffer_cinemri=temp_buffer_cinemri;
        min_rangex=rangex(1)*2;max_rangex=rangex(end)*2;
        rangex=min_rangex:max_rangex;
        min_rangey=rangey(1)*2;max_rangey=rangey(end)*2;
        rangey=min_rangey:max_rangey;
    elseif(upsample_factor==4)
        for i=1:size(buffer_cinemri,3)
            temp_buffer_cinemri(:,:,i)=interp2(buffer_cinemri(:,:,i),2,'cubic');
        end
        buffer_cinemri=temp_buffer_cinemri;
        min_rangex=rangex(1)*4;max_rangex=rangex(end)*4;
        rangex=min_rangex:max_rangex;
        min_rangey=rangey(1)*4;max_rangey=rangey(end)*4;
        rangey=min_rangey:max_rangey;
    end
end
%   buffer_cinemri=buffer_cinemri(:,:,ranget);
temp=strcat('Output/',int2str(slice),'/buffer_cinemri.mat');
save(temp,'buffer_cinemri');

cinemri=buffer_cinemri(:,:,:);   %%% DEV FOR TEST ON MARCH 18
 temp=strcat('Output/',int2str(slice),'/cinemri.mat');
 save(temp,'cinemri');
cinemri1=buffer_cinemri(:,:,:);   %%% DEV FOR TEST ON MARCH 18
%[cinemri_least_squares,shfts]=staticshft(cinemri1,[4 3 2 1], referenceFrame, numSkip,rangex,rangey,slice);

[RV LV]=FindLVRV(cinemri1);
x=floor(mean([RV(1),LV(1)]));
y=floor(mean([RV(2),LV(2)]));
% imagesc(cinemri1(:,:,13));colormap gray;axis image;
% [x,y]=ginput(1)
% x=floor(x);
% y=floor(y);
% close all
[shifts_out, cinemri_least_squares] = NeighboringTemporalSmoothingRegistration(cinemri1, x-40:x+40, y-40:y+40, 15, 15,0);
clear cinemri


% cinemri1=tmp_imgs(:,:,:);   %%% DEV FOR TEST ON MARCH 18
temp=strcat('Output/',int2str(slice),'/cinemri1.mat');
save(temp,'cinemri1');
temp=strcat('Output/',int2str(slice),'/cinemri_least_squares.mat');
save(temp,'cinemri_least_squares')
% done

% stage 3.1 - drawing RV region
% try
%     fname=strcat('Output/',int2str(slice),'/bw_rv.mat');
%     load(fname)
%     if(upsample_flag==1)
%         if(upsample_factor==2)
%             bw=interp2(bw,1,'cubic');
%         elseif(upsample_factor==4)
%             bw=interp2(bw,2,'cubic');
%         end
%     end
% catch    
    imagesc(cinemri1(:,:,referenceFrame)),colormap gray,brighten(0.3)
    disp('Draw RV region for input function')
    temp=strcat('Output/',int2str(slice),'/cinemri_least_squares.mat');
    load(temp)
       %[bw x y]=roipoly;
       close all
       [bw x]=auto_roi_mbr(cinemri_least_squares);
    temp=strcat('Output/',int2str(slice),'/bw_rv.mat');
save(temp,'bw','x','y');
% end
%done

g_sh=zeros(size(cinemri1,3),2,noi);

for MBR_iter_no=1:noi    
    tic
    
      temp=strcat('Output/',int2str(slice),'/cinemri_least_squares.mat');
      load(temp)
    cinemri1=cinemri_least_squares;
    % stage 3.2 - extracting curves
    disp(strcat('Doing: extracting curves - iteration #',int2str(MBR_iter_no)))    
    ii=find(bw);    
    bldcurve=zeros(size(cinemri1,3),1);    
    parfor i=1:size(cinemri1,3)
        bldcurve(i)=sum(sum(bw.*cinemri1(:,:,i)))/length(ii);
    end
    tmpImg=ones(size(cinemri1(:,:,1)));    
    [Y,X]=find(tmpImg);    
    nX = length(X);    
    tisscurves=zeros(nX,size(cinemri1,3));    
    parfor j = 1 : nX
        tisscurves(j,:) = cinemri1(Y(j), X(j), :);
    end    
    curves=[bldcurve';tisscurves];    
    temp=strcat('Output/',int2str(slice),'/curves.mat');
save(temp,'curves');
    % done
    
    % stage 3.4 - deltaSIcurves    
    disp(strcat('Doing: generating deltaSIcurves - iteration #',int2str(MBR_iter_no)))    
    si_curves = curves;
    [nRegs, nTimes]=size(curves);
    bldcurve=si_curves(1,:)';
    nRegs=nRegs-1;
    tisscurves = si_curves(2:nRegs+1,:);
    tisscurve=tisscurves';
    init_bld=sum(bldcurve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
    tiss_avgSpatial=0; counter=0;
    
    init_tiss=zeros(nRegs,1);
    
    parfor ireg=1:nRegs
        init_tiss(ireg)=sum(tisscurve(numSkip+1:numPreContrast_tiss+numSkip,ireg))/numPreContrast_tiss;
        if(init_tiss(ireg)>0)
            tiss_avgSpatial=tiss_avgSpatial+init_tiss(ireg);
            counter=counter+1;
        end
    end
    tiss_avgSpatial=tiss_avgSpatial/counter;    
    
    deltaSI_bldcurve= (bldcurve - init_bld);
    deltaSI_bldcurve=tiss_avgSpatial*deltaSI_bldcurve./(tiss_avgSpatial*ones(length(bldcurve),1)) ;
    temp=strcat('Output/',int2str(slice),'/init_bld.mat');
save(temp,'init_bld');
    
    deltaSI_tisscurve=zeros(nTimes,nRegs);
    parfor ireg=1:nRegs
        deltaSI_tisscurve(:,ireg)=tisscurve(:,ireg)-init_tiss(ireg);
        init_tiss(ireg);
        if(init_tiss(ireg)~=0)
            deltaSI_tisscurve(:,ireg)=tiss_avgSpatial*deltaSI_tisscurve(:,ireg)./(init_tiss(ireg)*ones(length(tisscurve(:,ireg)),1)) ;
        end
    end
    
    temp=strcat('Output/',int2str(slice),'/init_tiss.mat');
save(temp,'init_tiss');
    temp=strcat('Output/',int2str(slice),'/tiss_avgSpatial.mat');
save(temp,'tiss_avgSpatial');
    deltaSI_curves =[ deltaSI_bldcurve'; deltaSI_tisscurve'];  % so in same format as curves
    temp=strcat('Output/',int2str(slice),'/deltaSI_curves.mat');
save(temp,'deltaSI_curves');
    % done
    
    % 3.999 - Murase fitting
    disp(strcat('Doing: model fitting - iteration #',int2str(MBR_iter_no)))    
    murase_fitg(deltaSI_curves,slice);   %measured = fv*bld + ( conv( bld , k1*exp(-k2*(t-t0))) )    
    
    % Replacing with Murase fits
    disp(strcat('Doing: model image generation - iteration #',int2str(MBR_iter_no)))    
    cine_curv_fit(0,MBR_iter_no,Y,X,slice);
    % done
    
%      disp(strcat('Doing: least squares registration with model images - iteration #',int2str(MBR_iter_no)))
%      g_sh(:,:,MBR_iter_no)=min_lsg(rangex,rangey,slice,MBR_iter_no);    
    toc    
end
 temp=strcat('Output/',int2str(slice),'/g_sh.mat');
save(temp,'g_sh');



