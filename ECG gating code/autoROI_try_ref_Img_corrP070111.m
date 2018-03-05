

% 1. read in dicom like in 1.2, or
% read in cinemri, just recon with no header?

% find rv and lv centers.. 

% run self-gating

% won't be as great as would like - divide into systole and diastole, and 
% run non-rigid registration. 

% Check how that looks.  Then re-do recon, with only those frames,
% and with forward warping in regul. only.  (fixed, later try adaptive)



clear

saveflag=1   % change to 1 to write out peak locations for recon files
% consider:
% save this, change recon - reordering, might put all systole together,
%then all diastole. (if no contrast changes.., but even still might do 
%it.  
% And other option (likely not as good!) - recon, then modify curve from
% ROIs to just use peaks or valleys

raid1path='/v/raid1/ed/MRIdata/Cardiac/Verio/P070111/';
raid1pathGA='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P070111/';
islice=4;  % chnage also dir name below!!
islice=3

seriesName='CV_Radial7Off_flexICE_rcECG_pdTEST_reconILOT35(107000)'
%seriesName='SOS_centric_IL_tz_sparser0401_sat120_12deg_5ml(174000)'  

dirname=strcat(raid1pathGA,'ReconData/',seriesName,'/')

files=dir(strcat(dirname,'*.dcm'))
            for bb=1:length(files)
                tmpp=dicomread( strcat(dirname,files(bb).name));
                %tmpp=flipud(tmpp);
                %imgrr(:,:,bb)=rot90(tmpp);
                sos(:,:,bb)=double(tmpp);
            end
            

% choose frame # and ROI
% maybe 85, or 72.. 

% frameNum=72 %55 %72  %53
% figure; clf; imagesc(sos(:,:,frameNum))

fullcinemri1=sos(:,:,100:round(bb/2));

  [RV,LV] = FindLVRV(fullcinemri1,1);
  
  % just center on RV: 
  
  sx=size(fullcinemri1,1)
  sy=size(fullcinemri1,2)
  
%   temp = runningRV(injectionName);temp.point = [temp.point;RV];temp.index = [temp.index i]; runningRV(injectionName) = temp;
%             temp = runningLV(injectionName);temp.point = [temp.point;LV];temp.index = [temp.index i]; runningLV(injectionName) = temp; 
            windowWidth = sx/3.3;
            windowHeight = sy/3.3;
            rangex = round(RV(1) - windowWidth/2):round(RV(1) + windowWidth/2);
            rangey = round(RV(2) - windowHeight/2):round(RV(2) + windowHeight/2);
            rangey = round(RV(2) - windowHeight/2):round(RV(2) + windowHeight);
            
            sz_adjust=1.5
            rangex2 = round(RV(1) - sz_adjust*windowWidth/2):round(RV(1) + sz_adjust*windowWidth/2);
            rangey2 = round(RV(2) - sz_adjust*windowHeight/2):round(RV(2) + sz_adjust*windowHeight);
            RV(1) = max(RV(1),min(rangex)+1);
            RV(2) = max(RV(2),min(rangey)+1);
            
            
  text(RV(2)+10, RV(1)-15,'RV','Color',[1 1 1]);
            plot([min(rangey) max(rangey) max(rangey) min(rangey) min(rangey)],[min(rangex) min(rangex) max(rangex) max(rangex) min(rangex)]);
            RV = RV - [min(rangex) min(rangey)]; %#ok<NASGU>
            LV = LV - [min(rangex) min(rangey)]; %#ok<NASGU>
            
            
            
 frameNum=53
imagesc(sos(:,:,frameNum))
% 
% tmpname=strcat(raid1path,'ReconData/mat_files/refImg_polyCoords_slice',int2str(islice),'.frame',int2str(frameNum));
% if(~exist(tmpname,'file'))
%     disp('Use the mouse to trace the boundary contour of RV and LV, then press "Enter" or double click.')
%     set(gcf,'Renderer','zbuffer');   % to stop the blinking. Trying 3/22/05
%     [bw1,x1poly,y1poly]=roipoly;
%     polyCoords=[x1poly y1poly];
%     save(tmpname, 'polyCoords','-ascii');
% else
%     ss=load(tmpname);
%     x1poly=ss(:,1); y1poly=ss(:,2);
%     bw1=roipoly(sos(:,:,frameNum),x1poly,y1poly);
%     %keyboard
% end
% 
% hold
% %fill(x1poly,y1poly,'r');
% h1=line(x1poly,y1poly);
% set(h1,'Color',[0 1 0]); 
% 
% 
% keyboard

refImg=sos(rangex,rangey,frameNum);

% popd(refImg);
% keyboard

% 
% 
% 
% [I, J]=find(bw1);
% 
% bw1=double(bw1);
% 
% refImg=sos(:,:,frameNum).*bw1;    %(I,J,frameNum);
% 
% newcropI=min(I):max(I); 
% newcropJ=min(J):max(J);
% 
% refImg=refImg(newcropI,newcropJ);      
% 
% bwROI=bw1;
for ii=1:size(sos,3)
    %curve1(ii)=mean(mean( (sos(:,:,ii).*bwROI) ));
    curve1(ii)=mean(mean( (sos(rangex,rangey,ii)) ));
%    curveSYS(ii)=mean(mean( (selected_sos(:,:,ii).*bwROI) ));
end
%plot(curve)


% compute corr. as in Larson 2004, to give a gating signal:
numFrames=size(sos,3);

ref=refImg-mean(refImg(:));
fbar=0 %mean(sos,3);
% for ii=1:numFrames
%   sos(:,:,ii)=sos(:,:,ii)-fbar; 
% end

denom2=sum(sum(ref.^2));
figure(10); imagesc(ref)

%for ii=1:numFrames/2+60
for ii=1:numFrames 
    ii
    
    currentImg=sos(rangex,rangey,ii);
   % currentImg=currentImg(newcropI,newcropJ);
   
    dd(ii)=sum(sum(currentImg));
    bb(ii)=sum(sum(sos(:,:,ii)));
    %cc(ii)=xcorr2(currentImg,refImg)
    currentImg=currentImg - mean(currentImg(:));
    
    ee(ii)=sum(sum(norm(ref-currentImg)));
    
    num=sum(sum(currentImg.*ref));
    denom1=sum(sum(currentImg.^2));
    denom=sqrt(denom1.*denom2);
    corr_eq5(ii)=num/denom;
    
%     currentImg=sos(:,:,ii).*bw1;
%     currentImg=currentImg(newcropI,newcropJ);
    num=sum(sum(currentImg.*refImg));
    denom2=sum(sum(refImg.^2));
    denom=sqrt(denom1.*denom2);
    cc(ii)=num/denom;
    
%     figure(1)
%     imagesc(currentImg)
%     figure(2);
%     imagesc(currentImg-ref)
%     title(ii); colorbar
    %pause
    
    
end
    

 
dd_diff=diff(-dd);
peakloc=[];
for ii=1:length(dd_diff)-1
    if ( sign(dd_diff(ii))>sign(dd_diff(ii+1)) )
        peakloc=[peakloc ii+1];
    end
    
end
[ppvals, peaklocs2]=findpeaks(-dd, 'minpeakdistance',2)
selected_sos=sos(:,:,peaklocs2);
figNum=33;
figure(figNum); clf;
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
dd=dd./max(dd);
cc=cc./max(cc);
lw=1.6
plot(dd,'k','linewidth',lw); hold on
%plot(cc,'g','linewidth',lw)
plot(peaklocs2,dd(peaklocs2),'bo','linewidth',lw); hold on
%plot(locs,cc(locs),'mo','linewidth',lw); hold on
xlabel('Frame number')
legend('Sum in region ', 'Correlation in region ')
%whitebg('k')
%whitebg(figNum,[0.3 .3 0.3])


%title('dd, from sum ')
figure; 

imagesc(selected_sos(rangex,rangey,40))
figure;
popd(selected_sos(rangex2,rangey2,:),' -dd minpeak2 from sum roi ')

% -dd seems more systolic
peakloc=peaklocs2   % fewer frames, otherwise the same. 
fileout=strcat(raid1path,'ReconData/mat_files/sum_roi_sum_minusdd_sys_slice',int2str(islice),'.mat')
fileout=strcat(raid1path,'ReconData/mat_files/selectedFrames_sum_roi_minusdd_sys_slice',int2str(islice),'.mat')
if saveflag
    save(fileout, 'peakloc')
end



%adding 5/26/11, to write out the half closest to systole, then half
%closest to diastole.. 
systole_lowpeaks=peakloc;
[ppvals, peaklocsTmp]=findpeaks(dd, 'minpeakdistance',2)
diastole_highpeaks=peaklocsTmp;
new_diastole=[]; new_systole=[];
for ii=1:length(dd)
    % closest peak is:
           
    [minval_systole, indexsys]=min(abs(ii-systole_lowpeaks))
    [minval_diastole, indexdias]=min(abs(ii-diastole_highpeaks))
    
    if (minval_systole==0 || minval_diastole==0)
        if minval_systole==0
            new_systole=[new_systole ii];
        else
            new_diastole=[new_diastole ii];
        end
    else
        if ( abs(dd(ii)-dd(systole_lowpeaks(indexsys))) >= abs(dd(ii)-dd(diastole_highpeaks(indexdias))) )
            new_diastole=[new_diastole ii];
        else
            new_systole=[new_systole ii];
        end
    end
    
end

figure;
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
plot(dd)
hold on
plot(new_systole,dd(new_systole),'ko','linewidth',lw)
plot(new_diastole,dd(new_diastole),'ro','linewidth',lw)
xlabel('Frame number')
ylabel('Arbitrary units')
legend('Sum in region ', 'Closest to systole', 'Closest to diastole')
selected_sos=sos(:,:,new_systole);
figure; 
popd(selected_sos(rangex2,rangey2,:),' new systole from sum roi ')
selected_sos=sos(:,:,new_diastole);
figure; 
popd(selected_sos(rangex2,rangey2,:),' new diastole from sum roi ')
%fileout=strcat(raid1path,'ReconData/mat_files/roi_sum_systolic_ref53_recon82_slice',int2str(islice),'.mat')
%saveflag=1
fileout=strcat(raid1path,'ReconData/mat_files/',seriesName,'_roi_sum_dd_Nearsystolic_slice',int2str(islice),'.mat')
if saveflag
    peakloc=new_systole;
    save(fileout, 'peakloc')
end
fileout=strcat(raid1path,'ReconData/mat_files/',seriesName,'_roi_sum_dd_Neardiastolic_slice',int2str(islice),'.mat')
if saveflag
    peakloc=new_diastole;
    save(fileout, 'peakloc')
end
%saveflag=0

%keyboard


peakloc=peaklocs2   % fewer frames, otherwise the same. 
fileout=strcat(raid1path,'ReconData/mat_files/sum_roi_sum_minusdd_systolic2_slice',int2str(islice),'.mat')
if saveflag
    save(fileout, 'peakloc')
end




figure
 imagesc(squeeze (selected_sos(rangex,round(median(rangey)),:)))
 load /v/raid1/ed/src/matlab/Code_1.3_fromBrianDec10_2010/spect.cmap
 colormap(spect./max(spect(:)))
 set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
xlabel('Frame number after skipping frames')
 axis square  %image
 h=gca;
set(h,'YTickLabel',[])


[ppvals, peaklocs2]=findpeaks(dd, 'minpeakdistance',2)
selected_sos=sos(:,:,peaklocs2);
% figure; 
% popd(selected_sos(crops,crops,:),' +dd minpeak2 from sum roi ')

% make a function
figure; 
 %imagesc(squeeze (selected_sos(crops(1)+40:crops(1)+100,crops(1)+55,:)))
 imagesc(squeeze (sos(rangex,round(median(rangey))+12,130:250)))
 
 colormap(spect./max(spect(:)))
 axis square  %image
 set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
xlabel('Frame number')
title('+dd ')
% turn off y axis numbers:
h=gca;
set(h,'YTickLabel',[])


figure; imagesc(squeeze (selected_sos(rangex,round(median(rangey))+12,130:250)))
 
 colormap(spect./max(spect(:)))
 

% figNum=34;
% figure(figNum); clf;
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',14)
% dd=dd./max(dd);
% cc=cc./max(cc);
% plot(dd,'k'); hold on
% plot(cc,'g')
% plot(peaklocs2,dd(peaklocs2),'bo'); hold on
% plot(locs,cc(locs),'mo'); hold on
% xlabel('Frame number')
% legend('Sum in region ', 'Correlation in region ')

peakloc=peaklocs2;
fileout=strcat(raid1path,'ReconData/mat_files/sum_roi_plusdd_systolic2_slice',int2str(islice),'.mat')
if saveflag
    save(fileout, 'peakloc')
end





bb_diff=diff(bb);
bbpeakloc=[];
for ii=1:length(bb_diff)-1
    if ( sign(bb_diff(ii))>sign(bb_diff(ii+1)) )
        bbpeakloc=[bbpeakloc ii+1];
    end
    
end
bbselected_sos=sos(:,:,bbpeakloc);
figure; 
plot(bbpeakloc,bb(bbpeakloc),'ro'); hold on
plot(bb)
title('bb, from ALL sum ')
figure; 
popd(bbselected_sos(crops,crops,:),'bb not minpeak2 from sum all ')



% scale=256/double(max(selected_sos(:)))
% for ii=1:size(selected_sos,3)
%   ss(:,:,1,ii)=uint8(scale*selected_sos(:,:,ii));
% end
% figure; montage(ss(crop,crop,1,6:25),'Size',[2 10])


% scale=256/double(max(sos(:)))
% for ii=1:size(sos,3)
%   ss(:,:,1,ii)=uint8(scale*sos(:,:,ii));
% end
% figure; montage(ss(crop,crop,1,66:85),'Size',[2 10])


roiwallx=55:57
roiwally=78:79

tmpp=mean(sos(roiwallx,roiwally,:),2);
tmpp2=squeeze(mean(tmpp,1));
plot( squeeze(mean(tmpp,1)) );  title('roiwall, over time, likely not wall!')


% from /v/raid1/lchen/ecg  to get h
% plot(2.5*(0:1:1399),h(1:1400),'r','linewidth',1)
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',14)
% xlabel('msec')





%load P041411_timeAfter_280frames_slice1.mat
% this gives 280 point vector of timesAfter R-wave (pretty corrupt?)
% try looking at:
% skip first xx 
% nTs=32;
% nTf=171
% %nTf=210
% timeAfter=timeAfter(nTs:nTf);
% 
% % group less than 400 msec to start:
% firstfind_indices=find(timeAfter>(570-(islice-1)*45));
% tmp_timeAfter=timeAfter(firstfind_indices);
% 
% newIndex=find(tmp_timeAfter< (1160-(islice-1)*45))
% locs = firstfind_indices(newIndex)
% 
% selected_sos=sos(:,:,locs);
% figure; 
% popd(selected_sos,' from ECG file, all <400msec ')
% fileout=strcat(raid1path,'ReconData/mat_files/fromECG_300to580ms_40frames_slice',int2str(islice),'.mat')
% if saveflag
%     peakloc=locs;
%     save(fileout, 'peakloc')
% end

            