function dicomsorter_recursive(top_of_path2read,outputpath,seriesopt,sliceopt);
% DICOMSORTER is a function that sorts dicom data from a given input path
%   path2read based on their sequence/series of acquisition and slice
%   ordering. By default it only sorts on series number. Usage as follows:
%
%   dicomsorter(path2read), dicomsorter(path2read,1) - sorts on series number
%   dicomsorter(path2read, 1, 1) - sorts on series number and the
%   respective series into their slices
%
% See also SORTONSERIESNUMBER, SORTONSLICENUMBER
%
% Sathya Vijayakumar
% July 2009, UCAIR, University of Utah
% Made recursive and mkdir, EVRD 7/09

if nargin < 4
    sliceopt = 0;
end;
if nargin < 3 || (seriesopt == 0)
    seriesopt = 1;
end;
if nargin < 2 
    outputpath = pwd;
else
    if exist(outputpath)~=7
        mkdir(outputpath)
    end
end;
if nargin < 1
    errordlg('Sorry, cannot proceed without input directory');
    return;
end;
    
h = dir(top_of_path2read)

h2 = waitbar(0,strcat('Please wait, reading dicom headers for directory ', top_of_path2read) );

for ii=3:length(h)
   if (h(ii).isdir)
       dicomsorter_recursive( strcat(top_of_path2read,'/',h(ii).name), outputpath)
       d{ii-2}.SeriesNumber = -99;
   else

%for i = 3:length(h)
    try
        d{ii-2} = dicominfo(strcat(top_of_path2read,'/',h(ii).name));
    catch
        strcat(top_of_path2read,'/',h(ii).name,'  is not a dicom file, skip and continue ')  % seems a poor way to handle this
        % this likely only works if it is the first file listed in the
        % directory!
        d{ii-2}.SeriesNumber = -99;
    end
    
        
    waitbar(ii/length(h),h2);    
   end
end;
close(h2);

if (d{1}.SeriesNumber~=-99)
   sortonseriesnumber(d,top_of_path2read,outputpath);
end


