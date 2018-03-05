clear all
close all
clc

%filename='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P041012/RawData/CV_2D_RADIAL_ECG_10APR2012_10h57m.ecg'  %ecgout.ecg'

filename='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P061412/RawData/CV_2D_RADIAL_ECG_14JUN2012_04h38m.ecg'
%filename='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P021712/RawData/CV_2D_RADIAL_ECG_17FEB2012_10h19m.ecg'
%filename='/v/raid1/ed/MRIdata/Cardiac/Verio/P040111B/RawData/LY_2D_RADIAL1_01APR2011_08h11m.ecg'
%filename='/v/raid1/ed/MRIdata/Cardiac/Verio/P041411/RawData/LY_2D_RADIAL1_14APR2011_02h56m.ecg'
% odd, really short, prior to log. Like early dry runs?  Then nothing
% written out for actual stress study?? Or did I just not see and copy it?
%filename='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P030112/RawData/CV_2D_RADIAL_ECG_01MAR2012_04h20m.ecg'
%filename='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P030612/RawData/CV_2D_RADIAL_ECG_06MAR2012_08h46m.ecg'
%filename='/v/raid1/gadluru/MRIdata/Cardiac/Verio/P021712/RawData/CV_2D_RADIAL_ECG_17FEB2012_10h19m.ecg'


plotflag=1
[ecg, trigger, stats, comments]=read_ecg_data(filename, plotflag);

% is this from stats above:

LogStartMPCUTime= 38309290        
LogStopMPCUTime= 38392235    

% LogStartMPCUTime=61050085                  
% LogStopMPCUTime=61111477

% given start time, figure out # of samples to skip to sync ecg and start
% of acquisition. 

load P061412_MID91_test_timestamps_slice3_rest  % from running readMeasDataVB15_oneCoil_tmp.m
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

triggerTimesForStudy=find(trigger(numToSkip:lastSample))*2.5 + tt(1);  % will give vector of triggers, value of each is the time it trigger
% should start this a bit earlier, could trigger, then first image a second
% later... 
% 1 second earlier is 400 samples. Do say 800
%triggerTimesForStudy=find(trigger(numToSkip-800:lastSample))*2.5 + tt(1) - 800*2.5;

% first, find largest trigger time that is smaller than current image. 
% or, start with first trigger, find all images with time less than next
% trigger, calc times. 


% index_trigger=1;
% for ii=1:length(timestamps_OneperImage_oneslice)
%    while  timestamps_OneperImage_oneslice(ii) < triggerTimesForStudy(index_trigger)
% 	  index_trigger = index_trigger - 1;
%    end
% 	  
% 	timeFromTrigger(ii) = timestamps_OneperImage_oneslice(ii) - triggerTimesForStudy(index_trigger);
% 	  
% end
	
%tt=tt+halfRays*2.5;
timeAfter=zeros(1,length(timestamps_OneperImage_oneslice));
for ii=3:length(timestamps_OneperImage_oneslice)   %, find closest trigger that is a smaller time. Then find distance from that trigger. save image(imageNumber, triggerNumber, timeAfterTrigger)
   trigTime=find(triggerTimesForStudy<tt(ii))
   imgInfo{ii}.trigIndex=trigTime(end)
   imgInfo{ii}.timeAfterTrigger=tt(ii) - triggerTimesForStudy(trigTime(end));
   timeAfter(ii)=imgInfo{ii}.timeAfterTrigger;
end

%save P061412_timeAfter_slice3_rest timeAfter

% look at vcg
range=5000:10000
figure; plot(ecg(1,range),ecg(2,range),'-')
hold on; plot(ecg(1,trigger(range)),ecg(2,trigger(range)),'ro')

figure;
plot(ecg(1,numToSkip+7000:numToSkip+12000),ecg(2,numToSkip+7000:numToSkip+12000))
hold on
plot(ecg(1,trigger(numToSkip+7000:numToSkip+12000)),ecg(2,trigger(numToSkip+7000:numToSkip+12000)),'ro')

