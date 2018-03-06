function processGated_dev(raid1path, seriesName, first_series, numSlices,preBolusFlag, saveflag);


if preBolusFlag==1
    isub_range=[98]
else
    isub_range=[99]
end




 for islice=1:1:numSlices  % chnage also dir name below!!

iseries=first_series+(islice-1)*1000;

seriesName_tmp=[seriesName,num2str(iseries),')']

cd([raid1path,'/Processing/']);
numFrames=109;

call_sys_dias_modifyParFile(raid1path, islice, iseries, isub_range, numFrames);
 end
 
 
 
 
 
 function call_sys_dias_modifyParFile(raid1path, islice, iseries, isub_range, numFrames)

templateParFile=['series',int2str(iseries),'.slice',int2str(islice),'.par'] 
    
%peakloc_complement=setdiff(1:max(peakloc),peakloc)
raid1path=[raid1path, '/ReconData/mat_files/'];
if ~exist(raid1path, 'dir')
    mkdir(raid1path)
end

flag_found_framesToSelect=0;
    
for isub=isub_range  % low dose sys, low dose dias, high sys, high dias
      
   outfile=['series',int2str(iseries+isub),'.slice',int2str(islice),'.par']  % 91 is low dose diastole
   fidout=fopen(outfile,'w')
 
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
        tline
        %for each intervention check if this line has what we want
        %  if so write our own line
%        if(~isempty(strfind(tline,'framesToSelect')))
%             flag_found_framesToSelect=1;
 %           fprintf(fidout, 'framesToSelect=[7:',int2str(numFrames),'];\n');
            
%         if(~isempty(strfind(tline,'ranget')))
%             fprintf(fidout, 'ranget=[');
%             fprintf(fidout, '%d ', peakloc);
%             fprintf(fidout, ']\n');
            
%         elseif(~isempty(strfind(tline,'ranget')))
%             fprintf(fidout, 'ranget=[1:110]\n');   % a guess at low dose range
        if(~isempty(strfind(tline,'ranget')))
            if numFrames < 50         
                fprintf(fidout, ['ranget=[7:',int2str(numFrames),']\n']);
            else
                if isub==98        
                   fprintf(fidout, 'ranget=[7:40]\n');  
                else
                   fprintf(fidout, ['ranget=[41:',int2str(numFrames),']\n']);  
                end   
            end
        elseif(~isempty(strfind(tline,'scaleAIF')))
            if isub==98
                fprintf(fidout, 'scaleAIF=1\n');  
            else
                %fprintf(fidout, 'scaleAIF=10\n'); 
                fprintf(fidout, 'scaleAIF=1\n'); % going to stick with own AIF first for now 10/6/11
            end           
        elseif(~isempty(strfind(tline,'studyNum')))
            fprintf(fidout, 'studyNum=%d\n',iseries+isub);   % a guess at low dose range
        elseif(~isempty(strfind(tline,'seriesNumAIF')))
            if isub==98
                fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub);
            else
                %fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub-2); 
                fprintf(fidout, 'seriesNumAIF=%d\n', iseries+isub);  % going to stick with own AIF first for now 10/6/11
            end
        else
            fprintf(fidout, [tline '\n'] );
        end
    end
    
    fclose(fidin);
end

return
 