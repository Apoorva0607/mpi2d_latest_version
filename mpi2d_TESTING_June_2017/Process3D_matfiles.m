clear all
clc
close all
  infile="/home/pedgaonk/registration_3d/RegisteredData_P033117_MID73.mat";
load(infile);
%  a=new_reduced_non_k_space;
% infile='/v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P061416/ReconData/mat_files/3D_dual_slab/dia/PDs/Pre_Interp_3D_MID_58.mat';
% load(infile);
% b=new_reduced_non_k_space;
% for i=1:100
% if i<=13
% c(:,:,:,i)=b(:,:,:,i);
% else
% c(:,:,:,i)=a(:,:,:,i-13);
% end
% end


askedForParFileOverwrite=0;
%        image_data=double(abs(new_reduced_non_k_space).*1e7);

%                image_data=double(1e10.*abs(recon_data_inter));
%                image_data=double(abs(recon_data_inter));
Image=Data.SystoleAIF;
   image_data=double(1e10.*abs(Image));
if sum(size(Image))>470
    addon=0000;
    indexOfFoldersToDo = 1:size(image_data,3);
    st=size(image_data,4);
else
    addon=0000;
    indexOfFoldersToDo =1:size(image_data,4);
    st=size(image_data,3);
end


rangex=1:size(image_data,1);
rangey=1:size(image_data,2);
%     rangex=1:150;
%     rangey=1:100;

data_info=Data.Info;
save('data_info.mat','data_info');
load('data_info.mat');
TimeStamps_3D=data_info.TriggerTime';
save('timestamp_sys.mat','TimeStamps_3D')
%   time=timeStamp;
 load('timestamp_sys.mat');
%  time=TimeStamps_3D;
time = TimeStamps_3D-TimeStamps_3D(1);
%    time = timeStamp-timeStamp(1);
%
pd_end=data_info.NumberOfProtonDensityScans;
save('pd_end.mat','pd_end')
MyVal=0;


for i=indexOfFoldersToDo
    data.value = 0;
    if indexOfFoldersToDo<=2
        
       %%% For dual echo AIF where series 100 is 1st echo and 101 would be 2nd echo
        sliceNum=1;
        seriesNum=100*i+addon;
        header.SeriesDescription=100;
    
        
    else
        sliceNum=i;                  %%%%% For 3 D data
        seriesNum=1000*i+addon;
        header.SeriesDescription=1000;
    end


    

    %runningRV = []; runningLV = [];
    
    runningRV = containers.Map({''},{data}); remove(runningRV,'');
    runningLV = containers.Map({''},{data}); remove(runningLV,'');
    centerRV = containers.Map({''},{[]}); remove(centerRV,'');
    centerLV = containers.Map({''},{[]}); remove(centerLV,'');
    usePreviousRVLV = 0;
    MyVal = 0;   %really convoluted way of redoing the processing if there are a few slices that are unacceptable.  The redo is centered with the median LVRV position
    mymap = containers.Map({''},{data});
    remove(mymap,'');
    
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
    temp.infile = infile;
    temp.studyNum = seriesNum;
    temp.sliceNum = sliceNum;
    temp.seriesNumAIF = seriesNum;
    temp.sliceNumAIF = sliceNum;
    temp.ranget = ['1:' num2str(st)];
    temp.lastFrame = st;
    temp.framesToSelect=['1:' num2str(st)];  % added 9/20/11, need to check:
    temp.timeStampFile = ['timeStampSer' num2str(seriesNum) '.mat'];
   
    updateParFile(ParFileName,temp);
    
    %create the timestamp file
    timeStamp = time-time(1);
    save(temp.timeStampFile,'timeStamp');
    
    if sum(size(Image))>470
        cinemri1_tmp=squeeze(image_data(:,:,i,:));
    else
        %              cinemri1_tmp=imresize(image_data(:,:,:),[144, 144]);
        
        cinemri1_tmp=squeeze((image_data(:,:,:,i)));
        %               imagesc(cinemri1_tmp);
        %               mask=roipoly;
    end
    
    for count=1:st
        cinemri1_tmp1(:,:,count)=(((((cinemri1_tmp(:,:,count))))));
        
    end
    
    [RV,LV] = FindLVRV(cinemri1_tmp1,1);
    center_x=(RV(1)+LV(1))/2;
    center_y=(RV(2)+LV(2))/2;
    
    
    flag_pd=0;
    if ndims(image_data)==4 && flag_pd==0
        
%         cinemri1=cinemri1_tmp1(center_x-72:center_x+71,center_y-72:center_y+71,:);
        cinemri1=cinemri1_tmp1;
    else
        cinemri1=cinemri1_tmp1;
        
    end
    
    
    save(['Output/RVLVpos.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'RV','LV');
    save(['Output/cinemri1.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.mat'],'cinemri1')
    
    shfts=zeros(size(cinemri1,3),2);
    temp1=strcat(['Output/shiftsMAN.study' num2str(seriesNum) '.slice' num2str(sliceNum) '.txt']);
    mat2txt(temp1,shfts,1:st,1:st);
    %mymap(FoldersToDo{i}) = data;
end


disp('These are the guesses to where the heart is.  The RV point doesn''t matter, while the LV does for some stages.');
disp('If the cropping is incorrect please change the rangex and rangey in the par file');
disp('Here are the par files for your convienence');

