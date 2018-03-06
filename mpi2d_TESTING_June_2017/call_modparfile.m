function call_modparfile(parfile,param,value,addparam)  %DEV

%parfile--parfile to be edited
%param-parameter to be edited/added..string
%value-value of parameter..string
%addparam-flag to decide edit/add new

outfile=['tmpp_outt.par'];  
   fidout=fopen(outfile,'w');
 
   fidin = fopen(parfile);
    if(fidin < 0)
        disp('I could not find the template file. Create with stg=1.2. It should be here:');
        disp(parfile);
    end
    
if ~addparam
      
    while 1
        tline = fgetl(fidin);
        if ~ischar(tline)
            break;
        end
        
      if(~isempty(strfind(tline,param)))
            
                fprintf(fidout,'%s=%s\n',param,value); 
                %fprintf(fidout, '%framesToSelect=\n'); 
                    
                       
      else
            
            fprintf(fidout, [tline '\n'] );
      end
    end
 
    
else
    while 1
        tline = fgetl(fidin);
        
        if(~isempty(strfind(tline,param)))
             break
        end
               
        if ~ischar(tline)
            fprintf(fidout,'%s=%s\n',param,value);
            break
        end
        fprintf(fidout, [tline '\n'] );   
    end
    
    
end

    fclose(fidin);
    fclose(fidout);
    
    % 
    movefile('tmpp_outt.par', parfile);

return
