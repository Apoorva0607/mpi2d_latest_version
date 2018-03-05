function deltaSI_bldcurve=mpi_si2deltaSI(Mz0,curves,sliceNum,studyNum, outpath, numPreContrast_bld, numPreContrast_tiss, numSkip, scale_to_study1,ranget,flag_PD_end,pd_end)

% note: doesn't care if pixelwise or not, just read in a file of curves

if nargin==0,
    error('sliceNum argument needed for mpi_si2gdconc.');
end

if flag_PD_end==0
   si_curves = curves(:,pd_end+1+numSkip:end);
else
   si_curves = curves(:,1:end-pd_end);
end
   [nRegs, nTimes]=size(si_curves);
   bldcurve=si_curves(1,:)';
   nRegs=nRegs-1;
    tisscurves = si_curves(2:nRegs+1,:);
    tisscurve=tisscurves';
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEV
     
     init_bld=sum(bldcurve(1:numPreContrast_bld))/numPreContrast_bld;
%     
%      for i=1:1:size(tisscurves,1)
%    tisscurves(i,:)=(tisscurves(i,:)./(max(Mz0(i+1,:)))).*max(max(Mz0(2:end,:)));  %To correct for coil sensitivity PD SCALING   DEV
% %     
% %     %init_tiss(i)=sum(tisscurves_tmp(i,numSkip+1:numPreContrast_tiss+numSkip))./numPreContrast_tiss;
% %     init_tiss(i)=sum(tisscurves_tmp(i,1+numSkip:numPreContrast_tiss+numSkip))./numPreContrast_tiss;
% %     
% %     deltaSI_tisscurve(:,i)=tisscurves_tmp(i,:)-init_tiss(i);
% %         
%      end
       tiss_avgSpatial=0; counter=0;
    
       for ireg=1:nRegs
          init_tiss(ireg)=sum(tisscurves(ireg,1:numPreContrast_tiss))./numPreContrast_tiss;
          if(init_tiss(ireg)>0)
             tiss_avgSpatial=tiss_avgSpatial+init_tiss(ireg);
             counter=counter+1;
          end
       end
          tiss_avgSpatial=tiss_avgSpatial/counter;
    
    
        for ireg=1:nRegs   % first subtract off (doesn't matter each region different scales!) Then normalize for coils. Other way, need to normalize init_tiss...
            deltaSI_tisscurve(:,ireg)=tisscurves(ireg,:)-init_tiss(ireg);
                deltaSI_tisscurve(:,ireg)=tiss_avgSpatial*deltaSI_tisscurve(:,ireg)./(init_tiss(ireg)*ones(length(deltaSI_tisscurve(:,ireg)),1));
        end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEV
   
% subtract and divide for coil sens. (?) assuming bld same??
deltaSI_bldcurve= (bldcurve-init_bld);

% showcurves([deltaSI_bldcurve'; deltaSI_tisscurve'],'deltaSI scaled for
% coil'); %%%chamhe 09/08/17
% deltaSI_curves =[ deltaSI_bldcurve'; deltaSI_tisscurve'];  % so in same format as curves


return;



% adding from ../../3/code to make models!!!   9/12/04
flagUseModel=0;
%flagUseModel=1;
if flagUseModel~=1
 infilename=strcat(outpath,'cinemri1.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
  load(infilename);   % just loading to get size
[nrows ncols nTimes]=size(cinemri1);
tmpImg=zeros(nrows,ncols);
scale_img=ones(nrows,ncols);
init_img=zeros(nrows,ncols);

% next section just to calculate scale and init_images to write out
%nX     = length(X);
% should check nX==nRegs  (will for pixelwise, not for large regions!)
% for regions, need code that creates them...
flagPixelwise=0;
if flagPixelwise           %nX==nRegs
   disp('writing out scale and init images pixelwise')
   for i=1:length(X)
      scale_img(Y(i),X(i))=init_tiss(i)/tiss_avgSpatial;  % assume in same order
      init_img(Y(i),X(i))=init_tiss(i);
  end

else   % then not doing pixelwise
  disp('nada here')
  % dummy statement
end   % of if pixelwise
% adding 9/14/04 for write out scaleImg, just to use as dummy!
[X,Y,xCenter, yCenter, start_angle,bw1, bw2, bw3]= mpi_loadContours(tmpImg,sliceNum, studyNum, outpath);



% adding 5/12/04 just to do blasted 32 region case!!
% need to make this a function!!!!  Used in mpi_re_* and mpi_remake*
% stick in 5/12/04 for model with 32 regions:
 if (nRegs==32)
disp('doing 32 model')
   nSect=16;
   dAng  = 360 / nSect;
   endoFlag=1;
   [Xendo,Yendo,xCenter, yCenter, start_angle,bw1, bw3]= mpi_loadContours(tmpImg,sliceNum, studyNum,outpath,endoFlag);
   endoFlag=2;
   [Xepi,Yepi,xCenter, yCenter, start_angle,bw1, bw3]= mpi_loadContours(tmpImg,sliceNum, studyNum,outpath, endoFlag);

   for i = 1 : nSect
   nX=length(Xendo);
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle;
   Ang1 = i * dAng       + start_angle;
   Angles =  atan2(-(Yendo - yCenter), Xendo - xCenter) * 180/pi + 180;
   for j = 1 : nX
%      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)&(Angles(j) + 360 < Ang1)))
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)&(Angles(j) + 360 < Ang1)) | ((Angles(j)+720 > Ang0) & (Angles(j)+720 < Ang1)))

         cn = cn + 1;
         endoRoiX(cn) = Xendo(j);
         endoRoiY(cn) = Yendo(j);
      end

 end
   for n = 1 : cn
      scale_img(endoRoiY(n), endoRoiX(n)) = init_tiss(i)/tiss_avgSpatial;  % assume in same order
      init_img(endoRoiY(n), endoRoiX(n)) = init_tiss(i);
   end

   nX=length(Xepi);
   clear RoiX RoiY;
   cn = 0;
   Angles =  atan2(-(Yepi - yCenter), Xepi - xCenter) * 180/pi + 180;
   for j = 1 : nX
%      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)&(Angles(j) + 360 < Ang1)))
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)&(Angles(j) + 360 < Ang1)) | ((Angles(j)+720 > Ang0) & (Angles(j)+720 < Ang1)))
         cn = cn + 1;
         epiRoiX(cn) = Xepi(j);
         epiRoiY(cn) = Yepi(j);
      end
   end
   for n = 1 : cn
      scale_img(epiRoiY(n), epiRoiX(n)) = init_tiss(i+nSect)/tiss_avgSpatial;  % assume in same order
      init_img(epiRoiY(n), epiRoiX(n)) = init_tiss(i+nSect);
   end
  end
  end  % if 32




% do blood same way if pixelwise or not:
   [Ybld, Xbld] = find(bw1);
   for i=1:length(Xbld)
      scale_img(Ybld(i),Xbld(i))=1;
      init_img(Ybld(i),Xbld(i))=init_bld;
   end



outfilename=strcat(outpath,'scaleImg.study',int2str(studyNum),'.slice',int2str(sliceNum),'.prebld',int2str(numPreContrast_bld),'.pretiss',int2str(numPreContrast_tiss),'.float');
ff=fopen(outfilename,'w','ieee-be');
fwrite(ff,scale_img','float');
fclose(ff);
outfilename=strcat(outpath,'initImg.study',int2str(studyNum),'.slice',int2str(sliceNum),'.prebld',int2str(numPreContrast_bld),'.pretiss',int2str(numPreContrast_tiss),'.float');
ff=fopen(outfilename,'w','ieee-be');
fwrite(ff,init_img','float');
fclose(ff);

end  % if flagUseModel



showcurves([deltaSI_bldcurve'; deltaSI_tisscurve'],'deltaSI scaled for coil');
return;
