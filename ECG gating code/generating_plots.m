clear all
filename='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P061412/RawData/CV_2D_RADIAL_ECG_14JUN2012_04h57m.ecg'
plotflag=0
[ecg, trigger, stats, comments]=read_ecg_data(filename, plotflag);

firstfind_indices1=[3 14 16 19 21 24 27 29 32 34 36 38 41 43 46 48 50 53 55 57 60 62 65 67 70 72 75 77 80 82 85 88 90 93 95 98 101 103 106 109 111 114 116 119 122 125 127 130 132 135 137 140 142 145 148 150 153 155 158 161 163 166 169 171 174 176 179 182 184 187 189 192 194 197 200 202 205 207 210 213 215 217 220 222 225 228 231 234 236 239 242 245 248 251 254 257 ]
firstfind_indices2=[3 5 7 10 15 17 22 25 32 34 36 38 39 41 43 48 50 53 55 60 65 70 75 77 80 82 85 86 90 95 98 106 108 111 114 119 122 127 130 135 137 140 145 148 150 153 158 161 163 166 169 171 174 176 179 184 187 189 192 197 200 202 205 210 215 225 228 231 234 239 242 245 248 251 254 257 ];
%firstfind_indices1 are the frame numbers obtained after running self gating
%technique....Run autoProcess*.m... obtain the par files for systole and
%diastole and copy the frame numbers.
%firstfind_indices2 are the frame numbers obtained by manually runnibg the
%movie and noting down the frames for systole/diastole.
LogStartMPCUTime= 39430537        
LogStopMPCUTime= 39513505      

% LogStartMPCUTime=61050085                  
% LogStopMPCUTime=61111477

% given start time, figure out # of samples to skip to sync ecg and start
% of acquisition. 

load P061412_MID107_test_timestamps_slice3_stress  % from running readMeasDataVB15_oneCoil_tmp.m
%timestamps_perReadoutLine;
%halfRays=24/2
% can get TR from the readMeas timestamps
%TR=2.25

tt=2.5*timestamps_OneperImage_oneslice;  %- halfRays*2.5;  % now should match numbers in start time

% startAcqTime=timestamps_perReadoutLine(1)*2.5  - 12*2.5;   % if want midpoint 
% endAcqTime = timestamps_perReadoutLine(end)*2.5  - 12 *2.5;

numToSkip=(tt(1)-LogStartMPCUTime )/2.5

lastSample=numToSkip + (tt(end) - tt(1))/2.5 


ecgForstudy=ecg(1,numToSkip+1:lastSample);
triggersForstudy=trigger(numToSkip+1:lastSample);

totalTime_msec=(lastSample - numToSkip)*2.5


% for given slices, 24 readout rays per image, take say midpoint as the
% image time.  Can then look at how far after last trigger it is and bin
% it... 
% relate triggers to image number

triggerTimesForStudy=find(trigger(numToSkip:lastSample))*2.5 + tt(1);
timeAfter=zeros(1,length(timestamps_OneperImage_oneslice));
for ii=4:length(timestamps_OneperImage_oneslice)   %, find closest trigger that is a smaller time. Then find distance from that trigger. save image(imageNumber, triggerNumber, timeAfterTrigger)
   trigTime=find(triggerTimesForStudy<tt(ii));
   imgInfo{ii}.trigIndex=trigTime(end);
   imgInfo{ii}.timeAfterTrigger=tt(ii) - triggerTimesForStudy(trigTime(end));
   timeAfter(ii)=imgInfo{ii}.timeAfterTrigger;
end

firstfind_indices=find(timeAfter>(1) & timeAfter<120);

%plot(ecgForstudy)
hold on
%ecg_triggers=ecgForstudy.*triggersForstudy;
a=find(triggersForstudy);
%plot(a,ecgForstudy(a),'ro');
zoom xon
hold on
b=tt(firstfind_indices2);
b=(b-tt(1))/2.5;
b=round(b);
plot(b,ecgForstudy(b),'bx')
%plot(1:length(b),b,'go');
hold on
b1=tt(firstfind_indices);
b1=(b1-tt(1))/2.5;
b1=round(b1);
plot(b1,ecgForstudy(b1),'ro')
title('ECG Vs Manual');
figure
b1=tt(firstfind_indices2);
b1=(b1-tt(1))/2.5;
b1=round(b1);
plot(b1,ecgForstudy(b1),'ro')
hold on
b2=tt(firstfind_indices1);
b2=(b2-tt(1))/2.5;
b2=round(b2);
plot(b2,ecgForstudy(b2),'gx')
title('Self Vs Manual');

c=intersect(firstfind_indices,firstfind_indices1)

% a=timeAfter(firstfind_indices);
% hist(a,30)



