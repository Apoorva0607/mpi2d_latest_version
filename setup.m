function setup

raid1path='/v/raid1a/apedgaon/MRIdata/Cardiac/Prisma_version1/P082516_reg';
%  raid1pathDATA='/v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P102516/';
 raid1pathDATA1='/v/raid1a/ytian/MRIdata/Cardiac/3D_Data_For_Jason/P082516/processing';


  if ~isdir(raid1path)
      mkdir(raid1path)
  end
  cd(raid1path)
 
%   tmpstr=['ln -s ', raid1pathDATA, 'DicomData/* .']
%   unix(tmpstr)
%   cd ..
tmpstr=['ln -s ', raid1pathDATA1];
   unix(tmpstr)


%   tmpstr=['ln -s ', raid1pathDATA, 'RawData .'];
%    unix(tmpstr)

  mkdir('Processing_stress');
 
%  mkdir('ReconData')
%   cd('ReconData')
%     tmpstr=['ln -s ', raid1pathDATA, 'ReconData/* .'];
%     unix(tmpstr)
%     cd ..
    
%   !rm RawData
% %  mkdir('mat_files')
 cd('Processing_stress')
 mkdir('Output')
 cd('Output')
 cd ..
 cd ..
  mkdir('Processing_rest');
  cd('Processing_rest')
 mkdir('Output')
 cd('Output')

%  mkdir('mat_files')
 
