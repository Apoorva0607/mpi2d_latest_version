function processUngated(cinemri2,pd_frames,iseries, islice, preBolusFlag, saveflag)

if preBolusFlag==1
    isub_range=[90 91]
else
    isub_range=[92 93]
end

numFrames=call_autoROI(cinemri2,islice, saveflag,iseries);
% call_sys_dias_modifyParFile(cinemri2,pd_frames,islice, iseries, isub_range,numFrames);
return

function call_sys_dias_modifyParFile(cinemri2,pd_frames,islice, iseries, isub_range,numFrames)

templateParFile=['series',int2str(iseries),'.slice',int2str(islice),'.par']; 
    

flag_found_framesToSelect=0;
    
for isub=isub_range  
   if isub==90 || isub==92
       load(['Output/mat_files/listofselectedFrames_roi_sum_dd_Nearsystolic_series',int2str(iseries),'slice',int2str(islice),'.mat'])
   else
       load(['Output/mat_files/listofselectedFrames_roi_sum_dd_Neardiastolic_series',int2str(iseries),'slice',int2str(islice),'.mat'])
   end
   frames=peakloc(find(peakloc>5));
   cinemri1=cinemri2(:,:,frames);
   save(['Output/cinemri1.study' num2str(iseries+isub) '.slice' num2str(islice) '.mat'],'cinemri1');
   save(['Output/pd_frames.study' num2str(iseries+isub) '.slice' num2str(islice) '.mat'],'pd_frames');
   
   outfile=['series',int2str(iseries+isub),'.slice',int2str(islice),'.par'];
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
            fprintf(fidout, '%d ', peakloc);
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
        fprintf(fidout, '%d ', peakloc);
        fprintf(fidout, ']\n');
    end 
    fclose(fidin);
end

return



function numFrames=call_autoROI(cinemri2,islice,saveflag,iseries)

fullcinemri1=cinemri2(:,:,:);
soss=fullcinemri1;

 [LV,RV] = FindLVRV(fullcinemri1,1);  %interchanged LV RV DEV
  
  % just center on RV: 
  
  sx=size(fullcinemri1,1);
  sy=size(fullcinemri1,2);
  
            windowWidth = sx/2;
            windowHeight = sy/2;
            rangex = round(RV(1) - windowWidth/2):round(RV(1) + windowWidth/2);
            rangey = round(RV(2) - windowHeight/2):round(RV(2) + windowHeight);
            
            sz_adjust=1
            rangex2 = round(RV(1) - sz_adjust*windowWidth):round(RV(1) + sz_adjust*windowWidth/2);
            if round(RV(1) - sz_adjust*windowHeight) < 5
               rangex2=5:round(RV(1) + sz_adjust*windowHeight/2);
            end
            if round(RV(1) + sz_adjust*windowHeight/2) >sx
               rangex2=round(RV(1) - sz_adjust*windowWidth):sx-5;
            end
            
            rangey2 = round(RV(2) - sz_adjust*windowHeight/2)+15:round(RV(2) + sz_adjust*windowHeight+15);
            if 15+round(RV(2) - sz_adjust*windowHeight/2) < 5
               rangey2=5:round(RV(2) + sz_adjust*windowHeight+15);
            end
            if round(RV(2) + sz_adjust*windowHeight+15)>sy
               rangey2=round(RV(2) - sz_adjust*windowHeight/2)+15:sy-5;
            end
            RV(1) = max(RV(1),min(rangex)+1);
            RV(2) = max(RV(2),min(rangey)+1);
            
            
  text(RV(2)+10, RV(1)-15,'RV','Color',[1 1 1]);
  imagesc(fullcinemri1(:,:,50))
  hold on
  colormap gray
            plot([min(rangey2) max(rangey2) max(rangey2) min(rangey2) min(rangey2)],[min(rangex2) min(rangex2) max(rangex2) max(rangex2) min(rangex2)]);
            RV = RV - [min(rangex2) min(rangey2)]; %#ok<NASGU>
            LV = LV - [min(rangex2) min(rangey2)]; %#ok<NASGU>
  hold off          
         %keyboard   
          
numFrames=size(soss,3);
for ii=1:numFrames 
    currentImg=soss(rangex2,rangey2,ii);   %changed from rangex, EVRD 1/18/12, for P080411
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
selected_sos=soss(:,:,peaklocs2);
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
        abs((ii)-(systole_lowpeaks(indexsys)))
        if ( abs(dd(ii)-dd(systole_lowpeaks(indexsys))) <= (abs(dd(ii)-dd(diastole_highpeaks(indexdias)))))
            new_systole=[new_systole ii];
        else
            new_diastole=[new_diastole ii]; 
        end
    end
    
end

figure(700);
%set(gcf,'Color',[1 1 1])
%set(gca,'FontSize',14)
plot(dd)
hold on
lw=1.6;
plot(new_systole,dd(new_systole),'ko','linewidth',lw)
plot(new_diastole,dd(new_diastole),'ro','linewidth',lw)
%xlabel('Frame number')
%ylabel('Arbitrary units')
%legend('Sum in region ', 'Closest to systole', 'Closest to diastole')
selected_sos=soss(:,:,new_systole);
selected_sos=soss(:,:,new_diastole);

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


