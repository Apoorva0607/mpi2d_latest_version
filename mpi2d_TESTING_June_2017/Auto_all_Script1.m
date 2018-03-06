function Auto_all_Script1


list=dir(pwd);

% for j=5:6
%    cd(['/v/raid1a/dlikhite/MRIdata/Cardiac/Verio/Astellas_2013/For_SCMR/',list(j).name,'/Processing/']);
   folder=pwd;
  if(exist('mpi2d.par')) delete('mpi2d.par'); end
files = dir([folder '/*.par']);
for i=1:length(files)
   close all
   clc
    %copyfile(char(files(i).name),'mpi2d.par');
    eval(['!ln -s ',(files(i).name),' mpi2d.par'])
    disp(['i = ' num2str(i)]);

 try
% call_modparfile(char(files(i).name),'fixedDelay','1',0);        
% if i==1 || i==8
%     mpi2d stg=[3.2 3.4]
% else
%   mpi2d stg=[3.2 3.4 5.1]  
% end
                   mpi2d stg=[3.2]
%                           call_modparfile(char(files(i).name),'useDeltaSI','0',0);
%                            call_modparfile(char(files(i).name),'sliceNumAIF','1',0);
%                         call_modparfile(char(files(i).name),'seriesNumAIF','74099',0);
%                      call_modparfile(char(files(i).name),'numAzimuthalRegions','6',0);
%      %copyfile('mpi2d.par', char(files(i).name));
     !rm mpi2d.par
    
    catch
        !rm mpi2d.par
%           mail_error([pwd, files(i).name]);
% disp('ERROR')
% mail_error([list(j).name,' -s ',(files(i).name)])
        continue
    end
end 
   
  
%   
%  end
% mail  
return;




