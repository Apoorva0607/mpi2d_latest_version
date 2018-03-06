function cinemri1=func_getImagesFromDicomFolder(infolder);
% from mpi2d stg 0.2
    PathName=infolder
    DirectoryInfo = dir([PathName '/*']);
    if size(DirectoryInfo,1) == 0
        disp('there is nothing in this folder.  I am sorry');
        return;
    end
    %foreach HeaderInfo.SliceLocation
        %grab and parse times
        %sort times
    %pick lowest first time as first row in timeStamp
    %ascending from there
    SlicePositions = [];
    goodFilesFound = 0;
    for i=3:size(DirectoryInfo,1)   %first 2 entries are '.' and '..'
        try
            HeaderInfo = dicominfo([PathName '/' DirectoryInfo(i).name]);
            goodFilesFound = goodFilesFound+1;
            if isempty(find(SlicePositions == HeaderInfo.SliceLocation))
                SlicePositions(end+1) = HeaderInfo.SliceLocation;
            end
        catch
            
        end
    end
    
    SeriesNumber = zeros(goodFilesFound,1 );
    AcquisitionNumber = zeros(goodFilesFound,1 );
    InstanceNumber = zeros(goodFilesFound,1 );
    time = zeros(size(SlicePositions),floor(goodFilesFound./size(SlicePositions)) );
    lastTime = zeros(size(SlicePositions))+1;
    tframe = dicomread([PathName '/' DirectoryInfo(4).name]);  % changed to 4 EVRD  2/28/10
    cinemri1 = zeros(size(tframe,1),size(tframe,2),goodFilesFound );
    subset = -1;
    index2 = 1;
    for i=3:size(DirectoryInfo,1)   %first 2 entries are '.' and '..'
         try
%             if(strfind(char(DirectoryInfo(i).name),'subset') > 0)
%                 A = sscanf(char(DirectoryInfo(i).name),'subset%d');
%                 subset = A(1)
%             end
            DirectoryInfo(i).name
            %keyboard
            HeaderInfo = dicominfo([PathName '/' DirectoryInfo(i).name]);
            cinemri1(:,:,i-2) = dicomread([PathName '/' DirectoryInfo(i).name]);
        
            
            %cinemri1(i-1) = dicomread(HeaderInfo);
            index = find(SlicePositions == HeaderInfo.SliceLocation);
            SeriesNumber(i-2) = HeaderInfo.SeriesNumber;
            %parse them into seconds and bring zero them so start time is 0
            time(index,lastTime(index)) = str2num(HeaderInfo.AcquisitionTime);
            lastTime(index) = lastTime(index)+1;
            %timeStampSeries?.mat contains an nslic x frames matrix called timeStamp
            AcquisitionNumber(i-2) = HeaderInfo.AcquisitionNumber;
            InstanceNumber(i-2) = HeaderInfo.InstanceNumber;
        catch
            disp('note did not work for this file, likely not dicom')
        end
    end
end  % function
