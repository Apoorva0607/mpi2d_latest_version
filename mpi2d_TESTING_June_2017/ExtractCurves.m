function curves = ExtractCurves(cinemri1,sliceNum,studyNum,outpath,angle,numAzimuthalRegions,referenceFrame)
    %do a test on the angle to see if it's legacy(mod(angle,360)) or updated(mod(angle,2*pi))
%     if sliceNum==1
%         numAzimuthalRegions=2;
%     end
    
    legacy = 0;
    if(angle > 2*pi)
        legacy = 1;
        angle = mod((angle-90)*pi/180,2*pi);
    end
    
    [sx,sy,st] = size(cinemri1);
    %use the contours to recreate a mask as that's what the user would be
    %doing to recreate our results

%     outtermask = roipoly(cinemri1(:,:,tmpFrame),epi(:,1),epi(:,2));
%     innermask = roipoly(cinemri1(:,:,tmpFrame),endo(:,1),endo(:,2));
%     newMyo = (outtermask - innermask)>0;
%     newLV = roipoly(cinemri1(:,:,tmpFrame),blood(:,1),blood(:,2));
    %[r,c] = find(newLV>0);
    %newMeanLV = [mean(r) mean(c)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
figure(1)
tmpImg=mean(cinemri1(:,:,:),3); % not sure if this is needed 
tmpImg = tmpImg - min(tmpImg(:));
tmpImg = tmpImg / max(tmpImg(:));
RGBimg(:,:,1) = tmpImg;
RGBimg(:,:,2) = tmpImg;
RGBimg(:,:,3) = tmpImg;
imagesc(RGBimg);
axis('image')
    
endo_flag=0; % 0=Full 1=SubEndo 2=SubEpi   
    
[X,Y,xCenter, yCenter, start_angle,bw1, bw2, bw3, bw5]= mpi_loadContours(cinemri1,sliceNum, studyNum, outpath,endo_flag);  %DEV for midwall
if endo_flag==0
    newMyo = (double(bw2)-double(bw1))>0;
elseif endo_flag==1
    newMyo = (double(bw5)-double(bw1))>0;
elseif endo_flag==2
    newMyo = (double(bw2)-double(bw5))>0;
end
    


innermask=bw1;
outtermask=bw2;
newLV=bw3;
nX     = length(X);

dAng  = 360 / numAzimuthalRegions;
mask = zeros(sx, sy, numAzimuthalRegions);
Angles =  atan2(-(Y - yCenter), X - xCenter) * 180/pi + 180; 
% if worry about pixel sizes
%pixelSizeRatio=pixelsizeX/pixelsizeY; 
%pixelSizeRatio=1.8229/2.983
%pixelSizeRatio=1.0
mycolors = lines;

Angles =  atan2(-(Y - yCenter), (X - xCenter)) * 180/pi + 180; 
%%% NOTE - IF GET ERROR WITH RoiX, IS FROM START_ANGLE!!!  NOT ROBUST TO PICKING IN ANY  QUADRANT!!!
for i = 1 : numAzimuthalRegions
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle;
   Ang1 = i * dAng       + start_angle;
   for j = 1 : nX
%      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)))
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)) | ((Angles(j)+720 > Ang0) & (Angles(j)+720 < Ang1)))
         cn = cn + 1;
         RoiX(cn) = X(j); 
         RoiY(cn) = Y(j);
      end
   end
   
   for n = 1 : cn
      mask(RoiY(n), RoiX(n), i) = 1;
   end
   
   
%   for j = 1 : nFrames
%      t = img(RoiY, RoiX, j);
%      curve1(i,j) = mean(mean(t));
%   end
   
   figure(1);
   if(~exist('RoiX'))
       disp('Problem houston');
       keyboard
   end
   for j=1:length(RoiX)
       regionalColor = mycolors(i,:);
       myhsv = rgb2hsv(regionalColor);
       myhsv(3) = RGBimg(RoiY(j),RoiX(j),1);
       myrgb = hsv2rgb(myhsv);
       RGBimg(RoiY(j),RoiX(j),1) = myrgb(1);
       RGBimg(RoiY(j),RoiX(j),2) = myrgb(2);
       RGBimg(RoiY(j),RoiX(j),3) = myrgb(3);
   end
   imagesc(RGBimg), axis('image');
end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %use the myocardial mask to find the curves then the deltaSI curves
    [r,c] = find(outtermask>0);
    middle = [mean(r) mean(c)];
    [r,c] = find(newMyo>0);
    [theta,~] = cart2pol(r-middle(1), c - middle(2));
    theta = mod(theta - angle,2*pi);
    
    mycolors = reshape(lines(numAzimuthalRegions),1,numAzimuthalRegions,3);
    
    %subplot(2,2,3)   %EVRD
    %figure;
    myrgb = zeros(sx,sy,3);
    temp = mean(cinemri1(:,:,(referenceFrame-3):(referenceFrame+3)),3);
    temp = temp - min(temp(:));
    temp = temp/max(temp(:));
    myrgb(:,:,2) = temp;
    myrgb(:,:,3) = temp;
    myrgb(:,:,1) = temp;
    myregionalMask = zeros(sx,sy);
    labeled = floor(theta/(2*pi)*numAzimuthalRegions);
    %make it so that  region 3 and 4 are in the septal wall
    %(the region counter clockwise from the insertion angle is region 3, and the next clockwise region is region 4)
%     if(~legacy)
%         labeled = mod(labeled+2,6);
%     end
    for i=1:length(r)
        myregionalMask(r(i),c(i)) = labeled(i)+1;
    %    myrgb(r(i),c(i),:) = myrgb(r(i),c(i),:).*mycolors(1,labeled(i)+1,:);
    end
  %  imagesc(myrgb);
    
    for region=1:numAzimuthalRegions
        mymask = myregionalMask==region;
        [pointx,pointy] = find(mymask);
        pointx = mean(pointx);
        pointy = mean(pointy);
        text(pointy,pointx,num2str(region),'Color',[1 1 1]);
    end
    
    
    curves = zeros(st,numAzimuthalRegions);
    for region=1:numAzimuthalRegions
        regionalMask = myregionalMask == region;
        length(find(regionalMask==1))*0.23*0.23*0.8
        for t=1:st
            temp2 = cinemri1(:,:,t);
            a=temp2(regionalMask>0);
            curves(t,region) = mean(temp2(regionalMask>0));
        end
    end
    %other = load('/v/raid1/npack/Processing/P010710_nufft/Output/curves.study13.slice1.mat');
    %figure(2), subplot(2,2,1),plot(other.curves(2:end,:)'),subplot(2,2,2),plot(curves), subplot(2,2,3), plot(abs(other.curves(2:end,:)' - curves));
    %return;
    fifthHighestInAIF = zeros(st,1);
    for t=1:st
        temp = cinemri1(:,:,t);
        temp_pd=mean(cinemri1(:,:,1:4),3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        values1(:,:,t)=temp.*innermask;
        values_pd=temp_pd.*innermask;
        values2(:,:,t)=(values1(:,:,t)./values_pd).*(mean(mean(values_pd)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        meanInAIF(t) = mean(temp(newLV>0)); 
        
            
    end
max_aif=max(values2,[],3);
window_aif=find(max_aif>0.80*max(max(max_aif))&max_aif<=0.95*max(max(max_aif)));
temp_cine=cinemri1(:,:,referenceFrame);
temp_cine(window_aif)=0;
figure()
imagesc(temp_cine);
colormap gray
 for t=1:st
   temp=values1(:,:,t);
   aif_temp(:,t)=temp(window_aif);
 end
 
 %if sliceNum>1
 %aif=mean(aif_temp);
%  else
     aif=meanInAIF;
%  end
 
    figure()
    plot(aif,'r');
    hold on
    plot(meanInAIF)
    legend('autoAIF','drawnAIF')
     curves = vertcat(aif,curves');

    if(any(isnan(curves)))
        disp('Houston we have a curve problem');
        try
            if(exist('mpi2d.log'))
                fid = fopen('mpi2d.log','a'); 
            else
                fid = fopen('mpi2d.log','w'); 
            end
            fprintf(fid,'%s',['There were nans in the curves file typically implying that the contours led to invalid masks.  They were probably createdin stage 3.12./n']); 
            fclose(fid);
        catch ME %#ok<NASGU>

        end
        %keyboard;
    end
    %subplot(2,2,4), cla,hold on, plot(curves(1,:)','k');  %EVRD 12/15/10
    figure; plot(curves(1,:)','-*k');
    hold on
    for region=1:numAzimuthalRegions
        plot(curves(region+1,:)','-o','Color',squeeze(mycolors(1,region,:)));
    end
    title('Tissue+AIF curves');
    legend('AIF','1','2','3','4','5','6')
    
    return;
    
    