function makeParFileED(infolder, infile, outfolder, templateParFile, slice, subset, study,tmax)
    %load([infolder infile]);
    if(subset > 0)
        file = sprintf( 'rad.series%d.slice%d.subset%d.par',study,slice,subset);
    else
        file = sprintf('series%d.slice%d.par',study,slice);
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
        if(size(strfind(tline,'seriesNumAIF'))>0)
            fprintf(fid, 'seriesNumAIF=%d\n', study);
        elseif(size(strfind(tline,'sliceNumAIF'))>0)
            fprintf(fid, 'sliceNumAIF=%d\n', slice);
        elseif(size(strfind(tline,'studyNum'))>0)
            fprintf(fid, 'studyNum=%d\n', study);
        elseif(size(strfind(tline,'sliceNum'))>0)
            fprintf(fid, 'sliceNum=%d\n', slice);
        elseif(size(strfind(tline,'infile='))>0)
            fprintf(fid, 'infile=%s%s\n', infolder,infile);
        elseif(size(strfind(tline,'timeStampFile'))>0)
            fprintf(fid, 'timeStampFile=timeStampSer%d.mat\n',study );
        elseif(size(strfind(tline,'ranget'))>0)
            fprintf(fid, 'ranget=2:%d\n', tmax);
        elseif(size(strfind(tline,'lastFrame'))>0)
            fprintf(fid, 'lastFrame=%d\n', tmax-1);
        else 
            fprintf(fid, [tline '\n'] );
        end
    end
    fclose(fidin);fclose(fid);
end