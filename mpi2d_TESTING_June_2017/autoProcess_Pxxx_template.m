%function autoProcess_Pxxx_template(rangex,rangey,ranget,timeStampFile,infile,referenceFrame, numSkip,numPreContrast_bld,numPreContrast_tiss,seriesNum);
% to Process ungated perfusion.  9/19/11
% SET below parameters
% AFTER it is done, run 3.11 for all and fix, then mpi2d stg=[0.23 3.2 3.4 4.1
% 4.2]

% NOTE: Two functions put below in this file. Used to be separate scripts.


% reads in dicom like in 1.2, or
% read in cinemri, just recon with no header?
% find rv and lv centers.. 
% run self-gating



% SET THESE FOUR PARAMETERS

raid1path='/v/raid1/dlikhite/MRIdata/Cardiac/Verio/P030612/';
raid1pathDATA='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P030612/';
% first_seriesGated=
% seriesName='CV_Radial7Off_pd_ILOT35_goldenratio_Sept22params_full5.6ml('
% notew orth doing, just run 1.2 and mpi2d auto_all... 




% if ~isdir(raid1path)
%     mkdir(raid1path)
% end
% cd(raid1path)
% 
% %!ln -s /v/raid1/gadluru/MRIdata/Cardiac/Verio/P100411/ReconData/* .
% tmpstr=['ln -s ', raid1pathDATA, 'DicomData .']
% unix(tmpstr)
% 
% tmpstr=['ln -s ', raid1pathDATA, 'RawData .']
% unix(tmpstr)
% 
% mkdir('Processing')
% 
% mkdir('ReconData')
% cd('ReconData')
% tmpstr=['ln -s ', raid1pathDATA, 'ReconData/* .']
% unix(tmpstr)
% !rm mat_files
% mkdir('mat_files')
% cd('../Processing')


% run 1.2 for all files. 
% addpath('/v/raid1/ed/src/matlab/Code_1.5_Aug2011/')
%  mpi2d stg=1.2

 saveflag=1   % change to 1 to write out peak locations for recon files
% first_seriesUngated_lowdose=34000  % to 64000
  numSlices=5
% seriesNameUngated='CV_Radial7Off_flexICE_rcECG_pdTEST_reconILOT35('
% preBolusFlag=0
% processUngated(raid1path, seriesNameUngated, first_seriesUngated_lowdose, numSlices, preBolusFlag, saveflag);

first_seriesUngated=44000;
seriesNameUngated='CV_Radial7Off_golden_newSatprot_5.9ml_ungated('

% 
preBolusFlag=0
processUngated_dev(raid1path, seriesNameUngated, first_seriesUngated, numSlices, preBolusFlag, saveflag);

% run with own seriesNumAIF, 
% then re-run with some code to alter that series (and scale factor), and run stage 4.1 again. 
% Otherwise, need to pass in seriesNumAIF to processUngated.  which might be ok, then call it again, but it will 
% again find sys/dias frames etc..  too much.

% 
% cd([raid1path,'/Processing/']);
% select_SG_framesFlag=0;
% 
% first_seriesUngated_lowdose=61000 
% scaleAIF=10
% seriesNameUngated='CV_Radial7Off_pd_ILOT35_goldenratio_Sept22params_full5.6ml('
% 
% for islice=1:numSlices  % chnage also dir name below!!
% iseries=first_seriesUngated+(islice-1)*1000
% 
% %seriesName_tmp=[seriesNameUngated,num2str(iseries),')']
% 
% 
% sliceNumAIF=islice
% 
% %rFile(raid1path, islice, iseries, isub_range, seriesNumAIF, sliceNumAIF, scaleAIF);
% 
% % next part is from Auto_all_ForDualBolus
%     
% isub_range=[92 93]
%     for isub=isub_range  % low dose sys, low dose dias, high sys, high dias
%       %  if isub==90 || isub==92
%       parfile=['series',int2str(iseries+isub),'.slice',int2str(islice),'.par'] 
%       if(exist('mpi2d.par')) delete('mpi2d.par'); end
%      
%     % modify par file
%     seriesNumAIF=first_seriesUngated_lowdose+(islice-1)*1000+isub-2
%     call_modifyParFile(parfile, seriesNumAIF, sliceNumAIF, scaleAIF);
%     copyfile(parfile,'mpi2d.par');
%     
%       disp(['Processing ' parfile]);
%     
%       mpi2d stg=[3.4 4.1]
%     end
% 
% 
% end


    
    


