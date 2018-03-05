%ADDDICOMHEADER 
%   Function to add meta data from a given file, to a new dicom file that
%   is being created, with the specified tags modified.
%
%   Format: ADDDICOMHEADER(templateFilename, data_in,seriesNum, outname)
%       where, templateFilename is the name of the file whose dicom tags 
%       have to be copied
%       data_in is the input data that goes into the new dicom file thats
%       being written
%       seriesNum is the value of the series number tag that needs to be
%       replaced from the original dicom tag
%       seriesDesc is the series description that needs to be replaced from
%       the original dicom tag
%       instanceNum is the value of instance number that needs to be used
%       in the new dicom file
%       outname is the output file name of the new dicom file being
%       created.
function addDicomHeader(templateFilename, data_in, seriesNum, seriesDesc, instanceNum, outname)

%   Copyright UCAIR, University of Utah, 2007.
%   Author: Edward DiBella
%   Last Modified by: Sathya Vijayakumar Dec 6th 2007.


md=dicominfo(templateFilename);

md.SeriesDescription = seriesDesc;

md.SeriesNumber = seriesNum; 

data=data_in;
flatdata=reshape(data,1,prod(size(data)));


md.Height=size(data,1);  md.Width=size(data,2); 
outdata=double(data);

range=[];
if (isempty(range))
    range = [min(flatdata) max(flatdata)];
    scale = min(65535,max(flatdata));
end;
badData = find(~isfinite(outdata));

%scale= 65535/max(flatdata);
% scale= 511/max(flatdata);  % more like Siemens range, trying this 7/15/09 EVRD to see if quantization noise
%scale= max_int/max(flatdata);
scale=1;
outdata = uint16(scale*(outdata));

outdata(badData) = range(1);

md.BitDepth = 16;
md.BitsAllocated = 16;
md.BitsStored = 16;
   
md.SmallestImagePixelValue = scale*min(flatdata);
md.LargestImagePixelValue = scale*max(flatdata);
md.InstanceNumber = instanceNum;


%dicomwrite(squeeze(outdata),strcat(outname, num2strWithZeros(instanceNum),'.dcm') ,md);    % was outdata'
dicomwrite(squeeze(outdata),outname,md);    % was outdata'
    

function outnum=num2strWithZeros(num);
       if (num<=9)
          outnum=strcat('000',num2str(num));
       elseif ((9<num) && num<=99)
          outnum=strcat('00',num2str(num));
       elseif ((99<num) && num<=999)
          outnum=strcat('0',num2str(num));
       elseif ((999<num) && num<=9999)
          outnum=num2str(num);
       else
          disp('Problem with out of range number!!!')
       end
return

    
