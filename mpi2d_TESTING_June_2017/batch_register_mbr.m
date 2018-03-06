function test1
raid1path='/v/raid1/dlikhite/MRIdata/Cardiac/Verio/P012712/';
first_series=197000;
numSlices=10;
for islice=1:numSlices  % chnage also dir name below!!
iseries=first_series+(islice-1)*1000;
cd([raid1path,'/Processing/']);
parfile=['series',int2str(iseries),'.slice',int2str(islice),'.par'] 
      if(exist('mpi2d.par')) delete('mpi2d.par'); end
      copyfile(parfile,'mpi2d.par');
      mpi2d stg=[2.92]
end