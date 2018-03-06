clear all
load P061412_timeAfter_slice3_rest.mat      
%load the timeAfter*.mat file generated using the sync to ecg prog. 


%nTf=210
%timeAfter=timeAfter(nTs:nTf);
islice=1;
% group less than 400 msec to start:
firstfind_indices=find(timeAfter>(1) & timeAfter<170); %specify the range for the ecg window

tmp_timeAfter=timeAfter(firstfind_indices);
%firstfind_indices1=[2 4 7 9 11 13 15 18 22 24 27 30 32 34 38 42 45 47 49 52 54 56 58 60 63 67 69 71 74 77 79 81 83 85 87 89 91 97 99 101 103 105 108 112 115 117 120 122 125 128 130 133 135 137 142 144 146 148 150 152 154 156 159 162 166 170 172 174 176 178 180 182 186 189 191 196 199 202 204 208 210 212 214 218 220 223 225 228 230 232 236 238 240 244 248 ];
firstfind_indices2=[3 5 7 10 15 17 22 25 32 34 36 38 39 41 43 48 50 53 55 60 65 70 75 77 80 82 85 86 90 95 98 106 108 111 114 119 122 127 130 135 137 140 145 148 150 153 158 161 163 166 169 171 174 176 179 184 187 189 192 197 200 202 205 210 215 225 228 231 234 239 242];
flag=1
size(firstfind_indices)
load test_cine_P061412_s3_rest.mat   %load the cine file generated so as to test the gat
while flag==1;
for i=firstfind_indices;
    
    imagesc(fullcinemri1(89:203,82:196,i)); %the rangex and rangey are hardcoded using the .par file for that series and slice
    
    title(num2str(timeAfter(i)));
    colormap gray;
    pause
end

end



%To add dicom headers... comment out the aboe part and uncomment the lower
%part and specify the correct filenames.



% cinemri1=fullcinemri1(90:180,70:160,firstfind_indices);
% infile='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P082211/DicomData/CV_Radial7Off_pdTEST_ILOT35_ungated_STRESS(20)/'
% 
% 
% dicoms = dir([infile '*.dcm']);
%         if(isempty(dicoms))
%             disp('No dicoms found.  Cannot do Temporal Smoothing registration');
%             keyboard;
%             return;
%         end
%      for t=1:size(cinemri1,3)
%             
%         templatedirtxt='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P082211/DicomData/CV_Radial7Off_pdTEST_ILOT35_ungated_STRESS(20)/'
%         seriesnum=30000 + 44;
%         infilename=strcat('/v/raid1/gadluru/MRIdata/Cardiac/Verio/P082211/DicomData/CV_Radial7Off_pdTEST_ILOT35_ungated_STRESS(20)/',dicoms(firstfind_indices(t)).name)
%         outputdirtxt='/v/raid1/dlikhite/ECG gating code/P082211_self_gated_mod_stress_dia_s2';
%         outfilename=int2str(t);
%         indata=cinemri1(:,:,t);  
%         imagesc(indata)
%         outname = strcat(outputdirtxt,'/',outfilename);
%         seriesdesc='self gated mod stress dia s2'
%   
%         addDicomHeader(infilename, cinemri1(:,:,t), seriesnum,seriesdesc, t, outname);
%      end

