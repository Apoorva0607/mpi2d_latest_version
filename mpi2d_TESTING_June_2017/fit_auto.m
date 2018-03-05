function fit_auto


list=dir(pwd);

    folder=pwd;
  if(exist('mpi2d.par')) delete('mpi2d.par'); end
files = dir([folder '/*.par']);
for i=1:length(files)
   
  
   
    eval(['!ln -s ',(files(i).name),' mpi2d.par'])
    disp(['i = ' num2str(i)]);
    
 mpi2d stg=[4.2]
 
 !rm mpi2d.par
    
  
    end
end 