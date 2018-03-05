function file=makeParFile_forSubsets(infolder, outfolder, templateParFile, slice, subset, seriesNum, copyFromSeriesNum)
%modified from makeParFile to not modify anything in par file except infile (dicom folder point to)
% also use series # naming scheme
% EVRD 2/25/10

   
    %load([infolder infile]);
    if(subset > 0)
        file = sprintf( 'rad.series%d.slice%d.subset%d.par',seriesNum,slice,subset);
    else
        file = sprintf('series%d.slice%d.par',seriesNum,slice);
    end
    [fid message] = fopen([outfolder file] , 'w');
    if(fid < 0)
        folders = split(outfolder,'/');
        current = ['/' char(folders(1))];
        makeDirectories = 0;
        for j=2:size(folders,2)
            if(isempty(char(folders(j)))) 
                break;
            end
            if(exist(current,'dir') ~= 7)
                makeDirectories = 1;
            end
            if(makeDirectories > 0)
                mkdir(current,char(folders(j)));
                if(exist(current,'dir') ~= 7)
                    disp(['I tried to make a folder named:' char(folders(j)) 'but could not']);
                    disp(['Failed to open and process file:' infile]);
                    disp(['with the following message ' message]);
                    return;
                end
            end
            current = [current '/' char(folders(j))];
           
        end
        [fid message] = fopen([current '/' file] , 'w');
        if(fid < 0)
            disp('I tried to make the file but failed with this message:');
            disp(message);
            return;
        end
        %disp('I could not open the output folder.  Please make it for me.  Pretty pretty please');
    end
    %tmax = size( imgrr,3);
    fidin = fopen(templateParFile);
    if(fidin < 0)
        disp('I could not find the template file. It should be here:');
        disp(templateParFile);
    end
    while 1
        tline = fgetl(fidin);
        if ~ischar(tline)
            break;
        end
        %for each intervention check if this line has what we want
        %  if so write our own line
		if(size(strfind(tline,'seriesNum'))>0)
            fprintf(fid, 'seriesNum=%d\n', seriesNum);
        elseif(size(strfind(tline,'studyNum'))>0)
            fprintf(fid, 'studyNum=%d\n', seriesNum);
        elseif(size(strfind(tline,'sliceNum'))>0)
            fprintf(fid, 'sliceNum=%d\n', slice);
        elseif(size(strfind(tline,'infile='))>0)
            fprintf(fid, 'infile=%s\n', infolder);
        elseif(size(strfind(tline,'timeStampFile'))>0)
            fprintf(fid, 'timeStampFile=timeStampSer%d.mat\n',seriesNum );
		else 
            fprintf(fid, [tline '\n'] );
        end
    end
    
    fprintf(fid, 'copyFromSeriesNum=%d\n', copyFromSeriesNum);
            
            
            
    fclose(fidin);fclose(fid);
end