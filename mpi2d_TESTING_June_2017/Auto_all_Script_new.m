function Auto_all_Script_new


list=dir(pwd);

    folder=pwd;
  if(exist('mpi2d.par')) delete('mpi2d.par'); end
files = dir([folder '/*.par']);
for i=1:length(files)
   close all
   clc
   
    eval(['!ln -s ',(files(i).name),' mpi2d.par'])
    disp(['i = ' num2str(i)]);

%                                                          mpi2d stg=4.2
                                                     
                                                        mpi2d stg=[3.12 3.2 3.4 4.1]
                                                     
% % % % % % % %                                 
                  
                         
                     
     !rm mpi2d.par
    
  
    end
end 



