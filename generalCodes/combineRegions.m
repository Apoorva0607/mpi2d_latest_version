function [flowValuesOut] = combineRegions(flowValuesAHA)

% combineRegions  Take vector of flow values for apex to base in AHA seg.m model
% and combine to 3 territories
% of 18 regions, get LAD=   LCX=,   RCA=:w

% if 17 regions, combine:  


numRegionsOut=3;

% from Circ 2002 Cerqueira et al.
LADregions=[1,2,7,8,13, 14]; 
LCXregions=[5,6,11,12,16]; 
RCAregions=[3,4,9,10,15]; 

% modified for 18 regions (no cap and 6 regions in apical slice)
% LADregions=[1,2,7,8,13,14,15]; 
% LCXregions=[5,6,11,12]; 
% RCAregions=[3,4,9,10,16];     % note, not using 15

LADregions=[7,8]; 
LCXregions=[11,12]; 
RCAregions=[9,10];  


flowValuesOut(1)=mean(flowValuesAHA(LADregions));
flowValuesOut(2)=mean(flowValuesAHA(LCXregions));
flowValuesOut(3)=mean(flowValuesAHA(RCAregions));

return
