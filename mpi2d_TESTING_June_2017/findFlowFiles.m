function [ ktrans] = findFlowFiles(flowDir, study, slices,AIF_series,AIF_slice)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for ii=1:length(slices)
    
    
    % multiple AIF switch
    %if(length(AIFslice)==1)
        flowFile=[flowDir,'flowvalues.study',num2str((study(ii))),'.slice',num2str(slices(ii)),'.6.1_fixedDelay0.AIF_',num2str(AIF_series),'_',num2str(AIF_slice),'_1.txt.fermiFull.txt'];
        % 000 naming convention switch
%         if(~exist(flowFile,'file'))
%             flowFile=[flowDir,'flowvalues.study',num2str(study),'.slice',num2str(slices(ii)),'.6.1_fixedDelay0.AIF_',num2str(study),'_',num2str(AIFtype+AIFslice),'_1.txt.full.txt'];
%         end
%         if(~exist(flowFile,'file'))
%             flowFile=[flowDir,'flowvalues.study',num2str((study+slices(ii)-1)*1000),'.slice',num2str(slices(ii)),'.6.1_fixedDelay0.AIF_',num2str((study+AIFslice-1)*1000),'_',num2str(AIFtype+AIFslice),'_1.txt.full.txt'];
%         end
%         if(~exist(flowFile,'file'))
%             flowFile=[flowDir,'flowvalues.study',num2str(study*1000),'.slice',num2str(slices(ii)),'.6.1_fixedDelay0.AIF_',num2str((study+AIFslice-1)*1000),'_',num2str(AIFtype+AIFslice),'_1.txt.full.txt']
%         end
        
   % else
        %flowFile=[flowDir,'flowvalues.study',num2str((study+slices(ii)-1)*1000),'.slice',num2str(slices(ii)),'.6.1_fixedDelay0.AIF_',num2str((study+slices(ii)-1)*1000),'_',num2str(AIFtype+AIFslice(ii)),'_1.txt.full.txt'];
    %   flowFile=[flowDir,'flowvalues.study',num2str((study(ii))),'.slice',num2str(slices(ii)),'.6.1_fixedDelay0.AIF_',num2str((study(ii))*1000),'_',num2str(slices(ii)),'_1.txt.full.txt'];
        % 000 naming convention switch
%         if(~exist(flowFile,'file'))
%             flowFile=[flowDir,'flowvalues.study',num2str(study),'.slice',num2str(slices(ii)),'.6.1_fixedDelay0.AIF_',num2str(study),'_',num2str(AIFtype+AIFslice(ii)),'_1.txt.full.txt'];
%         end
%         if(~exist(flowFile,'file'))
%             flowFile=[flowDir,'flowvalues.study',num2str(study*1000),'.slice',num2str(slices(ii)),'.6.1_fixedDelay0.AIF_',num2str(study*1000),'_',num2str(AIFtype+AIFslice(ii)),'_1.txt.full.txt']
%         end
       
        
   % end
     if(~exist(flowFile,'file'))
            disp('have not found a filename that works for this yet, do by hand assign flowFile then return')
            %keyboard
%             ktrans=zeros(length(slices),6);
            break   % so ktrans will all be zero, use that as a flag to skip this person
     end
        
    ktrans_tmp=getKtrans(flowFile);
%     Vb(ii,:)=getVb(flowFile);
Vb=0;
    ktrans(ii,:)=ktrans_tmp(:);
end








end

