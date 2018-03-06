function [ecg, trigger, stats, comments]=read_ecg_data(plotflag)
% function [ecg, trigger, stats, comments]=read_ecg_data(filename, plotflag);
%
% function to read ecg data from Siemens PMU log

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

% default plot flag to be off


  filename = '';
  if length(filename) == 0
     [temp path] = uigetfile('*.ecg','Select File to Read');
     filename = [path temp(1:length(temp))];
  end



[p, ecg_filename]= fileparts(filename);

fid = fopen(filename);
i=0;
while ~feof(fid) 
   line = fgetl(fid);
   i=i+1;
   if isempty(line), break, end
   data{i}=line;
end
fclose(fid);
d1=data{1};

% comments between value 5002 (opening) and 6002 (closing)
index5002 = findstr(d1,'5002');
index6002 = findstr(d1,'6002');

% delete comments between value 5002 (opening) and 6002 (closing)
ecg=[];
for i=1:length(index6002)-1
    tmp= str2num(d1(index6002(i)+4:index5002(i+1)-1));
    len(i)=length(tmp);
    ecg=[ecg,tmp];
end
tmp=str2num(d1(index6002(end)+4:end));
len=cumsum([0,len]);
ecg=[ecg,tmp];

% parse out comments fields & tag the time index of message
tmp=[];
for i=1:length(index5002)
    comments{i}.index = len(i);
    comments{i}.message = d1(index5002(i)+4:index6002(i)-1);
%     tmp=strvcat(tmp,comments{i}.message);
end

PTABindex = [];
RELEARNindex = [];
for i=1:length(comments)
    switch comments{i}.message
        case {' MSGTYPE 110 ',' MSGTYPE 111 '} % relearn
            RELEARNindex = [RELEARNindex,comments{i}.index];
        case {' MSGTYPE 100 ',' MSGTYPE 101 ',' MSGTYPE 102 '} % table motion
            PTABindex = [PTABindex,comments{i}.index];
    end 
end       
RELEARNindex = round(RELEARNindex/2); % divide by 1 since ECG is dual channel
PTABindex = round(PTABindex/2);       % divide by 1 since ECG is dual channel

if nargout>=2
    stats=strvcat(data{2:end});
end


trigger=zeros(size(ecg));
start_triggers=find(ecg==5000);
for i=1:length(start_triggers)
    if start_triggers(i)==1
        trigger(start_triggers(i)+1)=1;
    else
        trigger(start_triggers(i)-1)=1;
    end
end

trigger=trigger(~or(ecg==5000,ecg==6000));

ecg=ecg(~or(ecg==5000,ecg==6000));


% check if 1st value is channel 1 or 2
if ecg(1) > 6000 % then must be channel 2
    ecg=ecg(2:end);
    trigger=trigger(2:end);
end
% check for even number of samples
if mod(length(ecg),2)~=0
    ecg=ecg(1:end-1);
    trigger=trigger(1:end-1);
end
ecg=reshape(ecg,[2 length(ecg)/2]);
trigger=reshape(trigger,[2 length(trigger)/2]);
trigger=or(trigger(1,:),trigger(2,:));

% remove bias
ecg(1,:)=ecg(1,:)-2048;
ecg(2,:)=ecg(2,:)-10240;

% delete last sample which is sometimes bad (??)
ecg=ecg(:,1:end-1);
trigger=trigger(:,1:end-1);


if plotflag>=1
    figure
    plot(ecg(1,:),ecg(2,:));
    hold on
    plot(ecg(1,trigger),ecg(2,trigger),'ro')
    axis equal
    t=[1:size(ecg,2)]*0.0025; % 400 samples/second
    ymax = 1.15 * max(abs(ecg(:)));
    h_fig=figure;
    h_axis(1)=subplot(2,1,1); plot(t,ecg(1,:));shg
    ecg_filename=strrep(ecg_filename,'_',' ');
    title(ecg_filename);
    hold on;plot(t(trigger),ecg(1,trigger),'ro');shg
    % delete values of RELEARNindex and PTABindex that are greater than ECG length
    RELEARNindex = RELEARNindex(RELEARNindex<length(ecg));
    PTABindex = PTABindex(PTABindex<length(ecg));
    if ~isempty(RELEARNindex)
        for i=1:length(RELEARNindex)
            plot([t(RELEARNindex(i)),t(RELEARNindex(i))],[-ymax ymax],'g')
        end
    end
    if ~isempty(PTABindex)
        plot(t(PTABindex),ecg(1,PTABindex),'ko')
    end
    axis([t(1) t(end) -ymax ymax]);
    h_axis(2)=subplot(2,1,2); plot(t,ecg(2,:));shg
    hold on;plot(t(trigger),ecg(2,trigger),'ro');shg
    if ~isempty(RELEARNindex)
        for i=1:length(RELEARNindex)
            plot([t(RELEARNindex(i)),t(RELEARNindex(i))],[-ymax ymax],'g')
        end
    end
    if ~isempty(PTABindex)
        plot(t(PTABindex),ecg(2,PTABindex),'ko')
    end
    axis([t(1) t(end) -ymax ymax]);
    xlabel('time (sec)')
    linkaxes(h_axis,'x');
    zoom xon
    pan(h_fig,'xon');
    
a=find(trigger);
b=mean(diff(a))*0.0025;
HR=(60/b)
end

%keyboard


return

% 
% 1. lPmuECGModePub                       // ECG control mode 
%         METHOD_NONE        = 0x01 
%         METHOD_TRIGGERING  = 0x02 
%         METHOD_GATING      = 0x04 
%         METHOD_RETROGATING = 0x08 
%         METHOD_SOPE        = 0x10 
%         METHOD_ALL         = 0x1E 
% 2. lPmuADPub                            // arrhythmia detection mode 
%         don't care 
% 3. iPmuHighPrioTriggerSignal            // high prio trigger signal 
%         SIGNAL_NONE        = 0x01, 
%         SIGNAL_EKG         = 0x02, 
%         SIGNAL_PULSE       = 0x04, 
%         SIGNAL_EXT         = 0x08, 
%         SIGNAL_RESPIRATION = 0x10, 
%         SIGNAL_ALL         = 0x1E, 
%         SIGNAL_EKG_AVF     = 0x20 
% 4. ulECGGateOnCountPub          // time stamp counter gate ON ECG 
%                                         // ON:  100 ms  / 2.5 = 40  
% 5. ulECGGateOffCountPub         // time stamp counter gate OFF ECG 
%                                         // OFF: 700 ms  / 2.5 = 280 
%                                  
% NOTE:
% If the latest Siemens VCG WIP (Improvements S073) has been installed on
% scanner, the outputs stats will include start and stop time for MDH and
% MPCU to allow better synchronization with the raw files
%
% LogStartMDHTime:                
% LogStopMDHTime:                 
% LogStartMPCUTime:               
% LogStopMPCUTime:   

% 
% below is a description of the header and footer data for the VB13A:
%  
% description of the addtional header data between 5002 and 6002:
%  
% eTriggerMethod: 1, minLimitCh1: 17, maxLimitCh1: 176, minLimitAVF: 3, maxLimitAVF: 38
%  
% eTriggerMethod: 1 means "standard VCG" trigger method, 2 means "channel I only", 3 means
%  "channel AVF only".
%  The remaining 4 values are just for internal debugging purposes, and can be ignored.
% description of the footer between 5003 and 6003:
% 
% short term statistics (current short term averages at the end of the data logging):
% 
% ECG Freq Per: 0 0
% PULS Freq Per: 0 0
% RESP Freq Per: 0 0
% EXT Freq Per: 0 0
% 
% long term statistics (long term averages etc. from the begin of the last protocol till the end of the data logging):
% 
% ECG Min Max Avg StdDiff: 0 0 0 0
% PULS Min Max Avg StdDiff: 610 2503 959 19
% RESP Min Max Avg StdDiff: 0 0 0 0
% EXT Min Max Avg StdDiff: 0 0 0 0
% NrTrig NrMP NrArr AcqWin: 49 0 0 921
% 


% additional info on comments fields

% // MSGTYPE in *.ecg file
% #define MSG_PTAB_STARTS_MOVING        100
% #define MSG_PTAB_STOPS_MOVING_INPOS   101
% #define MSG_PTAB_STOPS_MOVING_OUTPOS  102
% #define MSG_PTAB_POS_IS_IN            103
% #define MSG_PTAB_POS_IS_OUT           104
%  
% #define MSG_RELEARN_START             110
% #define MSG_RELEARN_STOP              111
% #define MSG_RELEARN_TIMEOUT           112
% #define MSG_RELEARN_ELECTRODEFAULT    113
%  
% #define MSG_LEARNED_DATA_WIP_FOPEN       210 // at fopen
% #define MSG_LEARNED_DATA_WIP_STOPSMOVING 211 // when PTAB stops
% #define MSG_LEARNED_DATA_WIP_RELEARN     212 // after relearn
%  
% #define MSG_LEARNED_DATA_VB_FOPEN       220 // at fopen
% #define MSG_LEARNED_DATA_VB_STOPSMOVING 221 // when PTAB stops
%  
%  
% I.e., 5002 MSGTYPE 100 6002 means, that at this time, the patient table (PTAB) started to move. 101 = PTAB stops moving in magnet, 102 = PTAB stops moving out of magnet. 103/104 PTAB position is in magnet / out of magnet - 103/104 are only used at the beginning of the logfile to have initial information on PTAB position. 110 and 111 indicate relearn start and relearn stop. 112 indicates, that relearning is stopped due to timeout (= no new waveform found within 20 sec = fall back to last waveforms). 113 marks, if relearing was stopped due to an electrode error when PTAB is out of the magnet. 
%  
% 210, 211, and 212 should have been followed(!?) by the learned data 
% - at logfile begin if PTAB is already inside magnet (210)
% - if PTAB stops inside the magnet and has been out when movement started (211)
% - after relearn (212)
%  
% Unfortunately, time was scarce in March last year so the comment code is written to the logfile but not the learned data.  ;-(  However, things are prepared. 
%  
%  
%  
% 220 and 221 indicate learned data of the VCG standard-algorithm
% - at logfile begin if is already inside magnet (220)
% - if PTAB stops inside the magnet and has been out when movement started (221)
%  
%  

                                        