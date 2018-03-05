function setup

% 


raid1pathDATA='/v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P092016/';
raid1path='/v/raid1a/apedgaon/MRIdata/Cardiac/Prisma/P092016/';
% first_seriesGated=
% seriesName='CV_Radial7Off_pd_ILOT35_goldenratio_Sept22params_full5.6ml('
% notew orth doing, just run 1.2 and mpi2d auto_all... 




  if ~isdir(raid1path)
      mkdir(raid1path)
  end
  cd(raid1path)
 

   mkdir('DicomData')
 cd('DicomData')
%   tmpstr=['ln -s ', raid1pathDATA, 'DicomData/* .']
%   unix(tmpstr)
  cd ..
%   tmpstr=['ln -s ', raid1pathDATA, 'RawData .']
%   unix(tmpstr)
 

  mkdir('Processing')
 
 mkdir('ReconData')
  cd('ReconData')
    tmpstr=['ln -s ', raid1pathDATA, 'ReconData/* .']
    unix(tmpstr)
 !rm mat_files
 mkdir('mat_files')
 cd('../Processing')
 mkdir('Output')
 cd('Output')
 mkdir('mat_files')
 
