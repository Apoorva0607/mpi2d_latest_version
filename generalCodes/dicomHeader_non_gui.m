% from Sathya's  function varargout = dicom_osirix_conv(varargin)


% set inpath, then only change 3 lines each time to do different series (if
% each a single slices as currently doing perfusion!!!  
inpath='/Volumes/HD2_Partition1/matlabResults/Radial/P070909/STCR_fromTaeHo'

templatedirtxt='/Volumes/Seagate931/Data/DICOM_FILES/P070909/CV_Radial7Off_3.9mladenoscan(26)'
seriesnum=10000 + 26
infilename='nufftP070909MID34_slice1_100iter_3coils_F1000_T900_S0'
outputdirtxt=templatedirtxt;
outfilename=infilename;
%--- maybe need to change below. depends on variable name in .mat file  -- %
load(strcat(inpath,'/',infilename) );   % gives imgrr  
% need a better way to figure out what the single variable name is!!
% maybe if exist im0 or sos or imgrr... 
indata=imgrr;   clear imgrr;
outname = strcat(outputdirtxt,'/',outfilename,'_');
seriesdesc=infilename(1:end-4);
func_dicomHeader(templatedirtxt, indata, outname, seriesnum, seriesdesc);


templatedirtxt='/Volumes/Seagate931/Data/DICOM_FILES/P070909/CV_Radial7Off_3.9mladenoscan(27)'
seriesnum=10000 + 27
infilename='nufftP070909MID34_slice2_100iter_3coils_F1000_T900_S0'
outputdirtxt=templatedirtxt;
outfilename=infilename;
load(strcat(inpath,'/',infilename) );   % gives imgrr  
indata=imgrr;   clear imgrr;
outname = strcat(outputdirtxt,'/',outfilename,'_');
seriesdesc=infilename(1:end-4);
func_dicomHeader(templatedirtxt, indata, outname, seriesnum, seriesdesc);


templatedirtxt='/Volumes/Seagate931/Data/DICOM_FILES/P070909/CV_Radial7Off_3.9mladenoscan(28)'
seriesnum=10000 + 28
infilename='nufftP070909MID34_slice3_100iter_3coils_F1000_T900_S0'
outputdirtxt=templatedirtxt;
outfilename=infilename;
load(strcat(inpath,'/',infilename) );   % gives imgrr  
indata=imgrr;   clear imgrr;
outname = strcat(outputdirtxt,'/',outfilename,'_');
seriesdesc=infilename(1:end-4);
func_dicomHeader(templatedirtxt, indata, outname, seriesnum, seriesdesc);





% now lexiscan:

templatedirtxt='/Volumes/Seagate931/Data/DICOM_FILES/P070909/CV_Radial7Off_3.9mladenoscan(54)'
outputdirtxt=templatedirtxt;
inpath='/Volumes/HD2_Partition1/matlabResults/Radial/P070909/STCR_fromTaeHo'
infilename='nufftP070909MID50_slice2_100iter_3coils_F1000_T900_S0'
outfilename=infilename;
seriesnum=10000 + 54

%--- maybe need to change below. depends on variable name in .mat file  -- %
load(strcat(inpath,'/',infilename) );   % gives imgrr  
% need a better way to figure out what the single variable name is!!
% maybe if exist im0 or sos or imgrr... 
indata=imgrr;
clear imgrr;
outname = strcat(outputdirtxt,'/',outfilename,'_');
seriesdesc=infilename(1:end-4);
func_dicomHeader(templatedirtxt, indata, outname, seriesnum, seriesdesc);



% 
% function func_dicomHeader(templatedirtxt, indata, outname, seriesnum, seriesdesc);
% % no help yet
% h = dir(templatedirtxt);
% slice_loc = [];
% n=0;
% for i = 3:42 %-- Here the assumption is for our scans, there are no more than 40 slices! -- %
%     h(i).name
%     
%     d = dicominfo(strcat(templatedirtxt,'/',h(i).name));
%     n = n+1;
%     if n == 1 || isempty(find(slice_loc == d.SliceLocation)) 
%         slice_loc = [slice_loc,d.SliceLocation];
%     else 
%         break;
%     end;
% end;
% 
% slice_loc2use = slice_loc(1);   % this is like 20.582
% 
% slice_match_indices = [];
% h1 = waitbar(0,'Please wait... Sorting through the template file directory');
% for i = 3:length(h)
%     waitbar(i/length(h),h1);
%     d = dicominfo(strcat(templatedirtxt,'/',h(i).name));
%     if d.SliceLocation == slice_loc2use
%         slice_match_indices = [slice_match_indices,i];
%     end;
% end;
% close(h1);
% 
% 
% % -- Start the dicom write process -- %
% if length(slice_match_indices) ~= size(indata,3)
%     %errordlg('Data size (3rd dimension) mismatch, please check input dataset! Will continue','Dataset size difference');
%     disp('Data size (3rd dimension) mismatch, please check input dataset! Will continue ');
% end;
% h2 = waitbar(0,'Please wait... Writing out the new dicom files');
% for k = 1:size(indata,3)
%     templateFilename = strcat(templatedirtxt,'/',h(slice_match_indices(k)).name);
%     
%     addDicomHeader(templateFilename, indata(:,:,k), seriesnum, ...
%         seriesdesc, k, outname);
%     waitbar(k/(size(indata,3)), h2);
% end;
% close(h2);
% return;

