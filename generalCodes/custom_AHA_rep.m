function [ threeFlows ] = custom_AHA_rep( ktrans)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% try to do AHA regions, then condense to 3 territories
% if 4 slices, apex to base, 
% first make base to apex,
% if gat==1
% ktrans.sys=ktrans.sys';
% else
% ktrans.dias=ktrans.dias';
% end
% if(size(ktrans.sys,2)==4)  % 4 slices 
%    tmp(:,4)=ktrans.sys(:,1);
%    tmp(:,3)=ktrans.sys(:,2);   % so will be base to apex
%    tmp(:,2)=ktrans.sys(:,3);
%    tmp(:,1)=ktrans.sys(:,4);
%  
ktrans=ktrans';
  tmp(:,1)=ktrans(:,1);
   tmp(:,2)=ktrans(:,2);   % so will be base to apex
   tmp(:,3)=ktrans(:,3);
 
   
   
   
   % try base, mid, then 4 apical, then 2 for cap. 
   %stressFlows=squeeze(reshape(tmp, [1,nRegs*nSlices])); 
   tmp=tmp(:);
   flows(1:12)=tmp(1:12);
   flows(13:16)=interp1(1:6, tmp(13:18),[1.5 3 4.5 6]);
   %flows.sys(17:18)=interp1(1:6, tmp(19:24),[2 4]);   
%    if (interpApicalSliceFlag==1)
%      stressFlows(13:16)=interp1(1:6, stressFlows(13:18),[1.5 3 4.5 6]);
%      stressFlows=stressFlows(1:16);
%    elseif (interpApicalSliceFlag==2)  % PET case, regions 17 and 18 are zero
%      stressFlows(13:18)=interp1(1:4, stressFlows(13:16),1:.6:4);
%    end
%    
% clear tmp
%    tmp(:,4)=ktrans.dias(:,1);
%    tmp(:,3)=ktrans.dias(:,2);   % so will be base to apex
%    tmp(:,2)=ktrans.dias(:,3);
%    tmp(:,1)=ktrans.dias(:,4);
%    
%    tmp=tmp(:);
%    flows.dias(1:12)=tmp(1:12);
%    flows.dias(13:16)=interp1(1:6, tmp(13:18),[1.5 3 4.5 6]);
%    flows.dias(17:18)=interp1(1:6, tmp(19:24),[2 4]);   
% 
% elseif(size(ktrans.sys,2)==5)
%    tmp(:,5)=ktrans.sys(:,1);
%    tmp(:,4)=ktrans.sys(:,2);
%    tmp(:,3)=ktrans.sys(:,3);   % so will be base to apex
%    tmp(:,2)=ktrans.sys(:,4);
%    tmp(:,1)=ktrans.sys(:,5);
%    % try base, mid, then 4 apical, then 2 for cap. 
%    %stressFlows=squeeze(reshape(tmp, [1,nRegs*nSlices])); 
%    tmp=tmp(:); 
%    flows.sys(1:6)=interp1(1:12,tmp(1:12),[1 3 5 7 9 11]);
%    flows.sys(7:12)=tmp(13:18);
%    flows.sys(13:16)=interp1(1:6, tmp(19:24),[1.5 3 4.5 6]);
%    flows.sys(17:18)=interp1(1:6, tmp(25:30),[2 4]);   
%    %stressFlows=stressFlows(1:18);
%    
%    clear tmp 
%    tmp(:,5)=ktrans.dias(:,1);
%    tmp(:,4)=ktrans.dias(:,2);
%    tmp(:,3)=ktrans.dias(:,3);   % so will be base to apex
%    tmp(:,2)=ktrans.dias(:,4);
%    tmp(:,1)=ktrans.dias(:,5);
%    % try base, mid, then 4 apical, then 2 for cap. 
%    %stressFlows=squeeze(reshape(tmp, [1,nRegs*nSlices])); 
%    tmp=tmp(:); 
%    flows.dias(1:6)=interp1(1:12,tmp(1:12),[1 3 5 7 9 11]);
%    flows.dias(7:12)=tmp(13:18);
%    flows.dias(13:16)=interp1(1:6, tmp(19:24),[1.5 3 4.5 6]);
%    flows.dias(17:18)=interp1(1:6, tmp(25:30),[2 4]);   
%    
%    
% else
%     disp('not implemented how to handle other than 4 or 5 slices ')
% end

    
threeFlows=combineRegions(flows)  % expects 18 regions AHA order
  
% threeFlows.dias=combineRegions(flows.dias)  % expects 18 regions AHA order

end

