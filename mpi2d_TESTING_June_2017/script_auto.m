function script_auto

list=dir(pwd);

    folder=pwd;

  if(exist('mpi2d.par')) delete('mpi2d.par'); end
files = dir([folder '/*.par']);
for i=1:length(files)
   close all
   clc
   
    eval(['!ln -s ',(files(i).name),' mpi2d.par'])
    disp(['i = ' num2str(i)]);


                                                       mpi2d stg=3.2
% % % % % % % % %                                 
%                   
                         call_modparfile(char(files(i).name),'useDeltaSI','0',0);
                         mpi2d stg=3.4
                         call_modparfile(char(files(i).name),'sliceNumAIF','1',0);
                        call_modparfile(char(files(i).name),'seriesNumAIF','100',0);
                         mpi2d stg=[ 4.1 5.1]
                     
     !rm mpi2d.par
    
  
    end
end 
   
  
  

                


