clear all
clc
close all

cd /v/raid11/gadluru/MRIdata/Cardiac/Prisma/P093016/ReconData/mat_files/2D/cSMS_perfusion/aif
list=dir('new_reduced_non_k_space_sms_MID_74_set_1.mat');
 load('/v/raid11/gadluru/MRIdata/Cardiac/Prisma/P093016/ReconData/mat_files/2D/cSMS_perfusion/aif/TimeStamps_MID_74.mat');



for num_files=1:length(list)

cd  /v/raid11/gadluru/MRIdata/Cardiac/Prisma/P093016/ReconData/mat_files/2D/cSMS_perfusion/aif
infile=list(num_files).name;    
load(infile);


cd /v/raid1a/dlikhite/MRIdata/Cardiac/Prisma/P093016/Processing_2D
askedForParFileOverwrite=0;



 
%image_data=double(abs(recon_data_inter)./1e5);


   addon=0;
   indexOfFoldersToDo = 3;
   
  
   time = (TimeStamps)/1e3;
  
   
   
   MyVal=0
    
 
for i=1:indexOfFoldersToDo
   
   if i==1
       image_data=double(abs(new_reduced_non_k_space_slice_1));
   elseif i==2
       image_data=double(abs(new_reduced_non_k_space_slice_2));
   else
       image_data=double(abs(new_reduced_non_k_space_slice_3));
   end
   
   rangex=1:size(image_data,1);
   rangey=1:size(image_data,2);
    st=size(image_data,3);  
      
    
    data.value = 0;
    sliceNum=i;
    
    seriesNum=100+i;
    header.SeriesDescription=1000*i;
    %runningRV = []; runningLV = [];
    
    runningRV = containers.Map({''},{data}); remove(runningRV,'');
    runningLV = containers.Map({''},{data}); remove(runningLV,'');
    centerRV = containers.Map({''},{[]}); remove(centerRV,'');
    centerLV = containers.Map({''},{[]}); remove(centerLV,'');
    usePreviousRVLV = 0;
    MyVal = 0;   %really convoluted way of redoing the processing if there are a few slices that are unacceptable.  The redo is centered with the median LVRV position
    mymap = containers.Map({''},{data});
    remove(mymap,'');

        

            %create the timestamp file
            

            if ndims(image_data)==4
            cinemri1_tmp0=squeeze(image_data(:,:,i,:));
            else
            cinemri1_tmp0=image_data;    
            end
            
            
            upsampleFactor=1;

            
            
            for count=1:st
                
               cinemri1_tmp(:,:,count)=interp2(cinemri1_tmp0(:,:,count),upsampleFactor-1,'cubic');
               cinemri1_tmp1(:,:,count)=(flipud((cinemri1_tmp(:,:,count))));
               
            end
            
            [RV,LV] = FindLVRV(cinemri1_tmp1,1);
             center_x=(RV(1)+LV(1))/2;
             center_y=(RV(2)+LV(2))/2;

            
            
            if ndims(image_data)==4

            cinemri1=cinemri1_tmp1(:,:,:);
            else
            cinemri1=cinemri1_tmp1(:,:,:);

            end
            
            
            save(['Output/RVLVpos.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'RV','LV');
            save(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'cinemri1')
            
            
            
     ParFileName = ['series' num2str(seriesNum) '.slice' num2str(sliceNum) '.par'];
            parFileMap(ParFileName) = header.SeriesDescription;
            if(exist(ParFileName) && ~askedForParFileOverwrite)
                choice = questdlg('I see some par files already here.  Should I use the rangex/y in them or overwrite with auto-cropping?','Par file usage','Reuse','Overwrite','Reuse');
                askedForParFileOverwrite = 1;
                switch choice
                    case 'Reuse'
                        ReadPar;
                    case 'Overwrite'
                end
            end
            if(~exist('Template.par'))
                disp(['I cannot find the Template.par file.  Please put it in a place I can find like (' myprocessingdir ') if you don''t have write permissions to the Code directory']);
                keyboard;
            end
            copyfile(which('Template.par'),ParFileName);
            clear temp;  
            
            temp.rangex = [num2str(min(rangex)) ':' num2str(max(rangex)-1)];  % Changed Dev for testing 6/13
            temp.rangey = [num2str(min(rangey)) ':' num2str(max(rangey)-1)];
            temp.infile = [infile];
            temp.studyNum = seriesNum;
            temp.sliceNum = sliceNum;
            temp.seriesNumAIF = seriesNum;
            temp.sliceNumAIF = sliceNum;
            temp.ranget = ['1:' num2str(st)];
            temp.lastFrame = st;
            temp.framesToSelect=['1:' num2str(st)];  % added 9/20/11, need to check: 
            temp.timeStampFile = ['timeStampSer' num2str(seriesNum) '.mat'];

            updateParFile(ParFileName,temp);        
            
            timeStamp = time-time(1);
            save(temp.timeStampFile,'timeStamp');
            
            
            
            
shfts=zeros(size(cinemri1,3),2);
temp1=strcat(['Output/shiftsMAN.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.txt']);
mat2txt(temp1,shfts,1:st,1:st);
            %mymap(FoldersToDo{i}) = data;
end


    disp('These are the guesses to where the heart is.  The RV point doesn''t matter, while the LV does for some stages.');
    disp('If the cropping is incorrect please change the rangex and rangey in the par file');
    disp('Here are the par files for your convienence');
   
end
