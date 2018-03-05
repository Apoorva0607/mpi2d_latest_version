% read in raw, recon with FBP and 3 subsets, and write out as dicom

% Note, views only shared within each 72 ray time frame (could extend to do 3 other time frames for 288 rays total). 
% Script needs dicom data unpacked by dicomsorter_recursive.m    
% 
% Takes in raw (.mat) k-space file and writes out recon with appropriate dicom headers, in 
% subdirectory to the Dicom data called Recon (can be altered. Note: don't write out to same Dicom directory,
% as subsequent uses of this program could get confused)

clear

scanner='Trio'
study='P010710'
MID=42
seriesnumStart=15
seriesnumEnd=16

MID=45
seriesnumStart=20
seriesnumEnd=21

MID=62
seriesnumStart=57
seriesnumEnd=58

% MID=63
% seriesnumStart=59
% seriesnumEnd=60


startFrame=1


weight_fidelity=1.0;
    weight_laplacian=0.9; % regularization parameter 1
    weight_TV=0.0; % regularization parameter 
    noi=100;
    stepsize=0.5;

    
    
% must also change filenames
filepath='/v/raid1/npack/MRIdata/Cardiac/Trio/P010710/ReconData/nufftRecons/combined/'




addpath /Users/ed/Documents/ShareWithWindows/matlab/generalCodes


% .02 rest study


for seriesnum=seriesnumStart:seriesnumEnd
  slice=seriesnum-seriesnumStart+1;
  filename=strcat('nufftP010710MID',int2str(MID),'_slice',int2str(slice),'_100iter_3coils_F1000_T900_S0.mat')
  load(strcat(filepath,filename))
%   template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',study,'/DicomData/CV_Radial7Off_triple_2.9ml(',int2str(seriesnum),')')
%   output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',study,'/ReconData/CV_Radial7Off_triple_2.9ml(',int2str(seriesnum*1000),')/');
  
% could find below by searching for one with seriesnum in the (), then
% parse it to give output directory... 
  template_dicom_directory=strcat('/v/raid1/bmatthew/MRIdata/Cardiac/',scanner,'/',study,'/DicomData/CV_Radial7Off_lexi_4.4ml(',int2str(seriesnum),')')
  output_dicom_directory=strcat('/v/raid1/ed/MRIdata/Cardiac/',scanner,'/',study,'/ReconData/CV_Radial7Off_lexi_4.4ml(',int2str(seriesnum*1000),')/')
  
  mkdir(output_dicom_directory)
  seriesdesc=strcat(study,'_series',int2str(seriesnum*1000),'_imgTCR_slice',int2str(slice))

  [sx sy sz]=size(imgrr);

            tmpg=imgrr((sx/4+1):(sx-sx/4),(sy/4+1):(sy-sy/4),:);clear imgrr
            imgrr=tmpg;

   outname=strcat(output_dicom_directory,'MID_',int2str(MID),'_');
   outname=strcat(outname,'TCR_slice_',int2str(slice),'_iter',int2str(noi),'_F',int2str(1000*weight_fidelity),'_T',int2str(1000*weight_laplacian),'_');
   
            

   func_dicomHeader(template_dicom_directory, imgrr, outname, seriesnum*1000, seriesdesc);
   clear imgrr
            %%% write out dicoms instead of saving .mats done ------ GA
            %%%
end


            