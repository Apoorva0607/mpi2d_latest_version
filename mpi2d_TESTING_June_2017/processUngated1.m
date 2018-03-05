function processUngated1(full_cinemri2,pd_frames,iseries, islice, preBolusFlag, saveflag,rangex_im,rangey_im,framesToSelect,ranget)

if preBolusFlag==1
    isub_range=[90 91]
else
    isub_range=[92 93]
end

numFrames=call_autoROI(full_cinemri2,islice, saveflag,iseries);
call_sys_dias_modifyParFile(full_cinemri2,pd_frames,islice, iseries, isub_range,numFrames,rangey_im,rangex_im,framesToSelect,ranget);

return

function call_sys_dias_modifyParFile(full_cinemri2,pd_frames,islice, iseries, isub_range,numFrames,rangey,rangex,framesToSelect,ranget)

templateParFile=['series',int2str(iseries),'.slice',int2str(islice),'.par']; 
    

flag_found_framesToSelect=0;
    
for isub=isub_range  
   if isub==90 || isub==92
       load(['Output/mat_files/listofselectedFrames_roi_sum_dd_Nearsystolic_series',int2str(iseries),'slice',int2str(islice),'.mat'])
   else
       load(['Output/mat_files/listofselectedFrames_roi_sum_dd_Neardiastolic_series',int2str(iseries),'slice',int2str(islice),'.mat'])
   end
   
   frames1=[1 2 3 4 peakloc];
   cinemri1_tmp=full_cinemri2(:,:,peakloc);
   
   close all
   if islice==1 
   cinemri2=cinemri1_tmp;
%   load(['Output/pd_frames.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat']);
   else
    run_registration(cinemri1_tmp,rangey,rangex,25, 1,5,5,iseries+isub);
    temp=strcat('Output/',int2str(iseries+isub),'/cinemri_curv_fit1.mat');
  load(temp);
  temp=strcat('Output/',int2str(iseries+isub),'/buffer_cinemri.mat');
 load(temp);
 cinemri2_tmp=ANTS_DMBR(buffer_cinemri,cinemri_curv_fit);
 
 run_registration(cinemri2_tmp,rangey,rangex,25, 1,5,5,iseries+isub+10);
    temp=strcat('Output/',int2str(iseries+isub+10),'/cinemri_curv_fit1.mat');
  load(temp);
  temp=strcat('Output/',int2str(iseries+isub),'/buffer_cinemri.mat');
 load(temp);
 cinemri2=ANTS_DMBR1(buffer_cinemri,cinemri_curv_fit);
 
 
 
 !rm -f *mhd
 !rm -f *raw
 !rm -f *gz
 !rm -f *txt
   end
 cinemri1_nup=zeros(size(cinemri2,1),size(cinemri2,2),size(cinemri2,3)+4);
 cinemri1_nup(:,:,1:4)=pd_frames(:,:,:);
 cinemri1_nup(:,:,5:end)=cinemri2;
 %save(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'cinemri1');
% load(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat']);
upsampleFactor=2;
for i=1:1:size(cinemri1_nup,3)
cinemri1(:,:,i)=interp2(cinemri1_nup(:,:,i),upsampleFactor-1,'cubic');
end


   save(['Output/cinemri1.study' num2str(iseries+isub) '.slice' num2str(islice) '.mat'],'cinemri1');
   shfts=zeros(size(cinemri1,3),2);
temp=strcat(['Output/shiftsMAN.study' num2str(iseries+isub) '.slice' num2str(islice) '.txt']);
mat2txt(temp,shfts,framesToSelect,ranget);
%    save(['Output/pd_frames.study' num2str(iseries+isub) '.slice' num2str(islice) '.mat'],'pd_frames');
   
   outfile=['Output/series',int2str(iseries+isub),'.slice',int2str(islice),'.par'];
   fidout=fopen(outfile,'w');
 
   fidin = fopen(templateParFile);
    if(fidin < 0)
        disp('I could not find the template file. Create with stg=1.2. It should be here:');
        disp(templateParFile);
    end
    while 1
        tline = fgetl(fidin);
        if ~ischar(tline)
            break;
        end
           
        if(~isempty(strfind(tline,'framesToSelect')))
            flag_found_framesToSelect=1;
            fprintf(fidout, 'framesToSelect=[');
            fprintf(fidout, '%d ', frames1);
            fprintf(fidout, ']\n');
        elseif(~isempty(strfind(tline,'ranget')))
                fprintf(fidout, ['ranget=[1:',int2str(numFrames),']\n']);
        elseif(~isempty(strfind(tline,'scaleAIF')))
            if isub==90 || isub==91
                fprintf(fidout, 'scaleAIF=1\n');  
            else
                fprintf(fidout, 'scaleAIF=1\n'); % going to stick with own AIF first for now 10/6/11
            end           
        elseif(~isempty(strfind(tline,'studyNum')))
            fprintf(fidout, 'studyNum=%d\n',iseries+isub);   % a guess at low dose range
        elseif(~isempty(strfind(tline,'seriesNumAIF')))
            if isub==90 || isub==91
                fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub);
            else
                
                fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub);  % going to stick with own AIF first for now 10/6/11
            end
        else
            fprintf(fidout, [tline '\n'] );
        end
    end
    if flag_found_framesToSelect==0
        fprintf(fidout, 'framesToSelect=[');
        fprintf(fidout, '%d ', frames1);
        fprintf(fidout, ']\n');
    end 
    fclose(fidin);
end

return



function numFrames=call_autoROI(full_cinemri2,islice,saveflag,iseries)

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
%     else
%         abs((ii)-(systole_lowpeaks(indexsys)))    %COMMENT FOR ONLY PEAKS
%         if ( abs(dd(ii)-dd(systole_lowpeaks(indexsys))) <= (abs(dd(ii)-dd(diastole_highpeaks(indexdias)))))
%             new_systole=[new_systole ii];
%         else
%             new_diastole=[new_diastole ii]; 
%         end
    end
    
end

figure(700);
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
plot(dd)
hold on
lw=1.6;
 new_systole=new_systole(find(new_systole>=5));
 new_diastole=new_diastole(find(new_diastole>=5));

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

return


