function gating_standalone(full_cinemri2)

islice=1;
iseries=10000;
saveflag=1;

fullcinemri1=full_cinemri2(:,:,:);
soss=fullcinemri1;

 [RV,LV] = FindLVRV(fullcinemri1,1);  

  % just center on RV: 
  
   sy=size(fullcinemri1,1);
   sx=size(fullcinemri1,2);
%   RV=[round(sx/2),round(sy/2)];
%  LV=[round(sx/2),round(sy/2)];
             windowWidth = sx/1.8;
             windowHeight = sy/2;
             H_cen(2)=(RV(2)+LV(2))/2;
             H_cen(1)=(RV(1)+LV(1))/2;
             rangex = round(H_cen(2) - windowWidth/2):round(H_cen(2) + windowWidth/2);
             rangey = round(H_cen(1) - windowHeight/2):round(H_cen(1) + windowHeight/2);
% rangex=30:sx-30;
% rangey=30:sy-30;
% 

rangex2=rangex;
rangey2=rangey;
            sz_adjust=1
%             rangex2 = round(RV(1) - sz_adjust*windowWidth):round(RV(1) + sz_adjust*windowWidth/2);
            if round(H_cen(2) - sz_adjust*windowWidth/2) < 5
               rangex2=5:round(H_cen(2) + sz_adjust*windowHeight/2);
            end
            if round(H_cen(2) + sz_adjust*windowWidth/2) >sx
               rangex2=round(H_cen(2) - sz_adjust*windowWidth/2):sx-5;
            end
            
%            rangey2 = round(RV(2) - sz_adjust*windowHeight/2):round(RV(2) + sz_adjust*windowHeight);
            if round(H_cen(1) - sz_adjust*windowHeight/2) < 5
               rangey2=5:round(H_cen(1) + sz_adjust*windowHeight/2);
            end
            if round(H_cen(1) + sz_adjust*windowHeight/2)>sy
               rangey2=round(H_cen(1) - sz_adjust*windowHeight/2):sy-5;
            end

rangex=rangex2;
rangey=rangey2;            
            
   imagesc(fullcinemri1(:,:,50))
   hold on
   colormap gray
   plot([min(rangex) max(rangex) max(rangex) min(rangex) min(rangex)],[min(rangey) min(rangey) max(rangey) max(rangey) min(rangey)]);

          
numFrames=size(soss,3);
for ii=1:numFrames 
    currentImg=soss(rangey,rangex,ii);

    dd(ii)=sum(sum(currentImg));
 
end
  

 
dd_diff=diff(-dd);
peakloc=[];
for ii=1:length(dd_diff)-1
    if ( sign(dd_diff(ii))>sign(dd_diff(ii+1)) )
        peakloc=[peakloc ii+1];
    end
    
end
[ppvals, peaklocs2]=findpeaks(-dd, 'minpeakdistance',1);

figNum=33;
figure(figNum); clf;
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
dd=dd./max(dd);

peakloc=peaklocs2;   % fewer frames, otherwise the same. 
fileout=strcat('Output/mat_files/listofselectedFrames_sum_roi_minusdd_sys_series',int2str(iseries),'slice',int2str(islice),'.mat');

if saveflag
    save(fileout, 'peakloc')
end


systole_lowpeaks=peakloc;
[ppvals, peaklocsTmp]=findpeaks(dd, 'minpeakdistance',1);
diastole_highpeaks=peaklocsTmp;
new_diastole=[]; new_systole=[];
for ii=1:length(dd)
           
    [minval_systole, indexsys]=min(abs(ii-systole_lowpeaks));
    [minval_diastole, indexdias]=min(abs(ii-diastole_highpeaks));
    
    if (minval_systole==0 || minval_diastole==0)
        if minval_systole==0
            new_systole=[new_systole ii];
        else
            new_diastole=[new_diastole ii];
        end
    else
        abs((ii)-(systole_lowpeaks(indexsys)))    %COMMENT FOR ONLY PEAKS
        if ( abs(dd(ii)-dd(systole_lowpeaks(indexsys))) <= (abs(dd(ii)-dd(diastole_highpeaks(indexdias)))))
            new_systole=[new_systole ii];
        else
            new_diastole=[new_diastole ii]; 
        end
     end
    
end

figure(700);
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
plot(dd)
hold on
lw=1.6;
%  new_systole=new_systole(find(new_systole>=10));
%  new_diastole=new_diastole(find(new_diastole>=10));

plot(new_systole,dd(new_systole),'ko','linewidth',lw)
plot(new_diastole,dd(new_diastole),'ro','linewidth',lw)
%xlabel('Frame number')
%ylabel('Arbitrary units')
%legend('Sum in region ', 'Closest to systole', 'Closest to diastole')


fileout=strcat('Output/mat_files/listofselectedFrames_roi_sum_dd_Nearsystolic_series',int2str(iseries),'slice',int2str(islice),'.mat')
if saveflag
    peakloc=new_systole;
    save(fileout, 'peakloc')
end
fileout=strcat('Output/mat_files/listofselectedFrames_roi_sum_dd_Neardiastolic_series',int2str(iseries),'slice',int2str(islice),'.mat')
if saveflag
    peakloc=new_diastole;
    save(fileout, 'peakloc')
end