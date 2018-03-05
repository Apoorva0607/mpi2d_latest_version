function [ gdcurves T1  ] = conc_gd1(curves,Mz0,sliceNum,seriesNum,pd_end,data_info)
flag_pd=1; %%%% flag_pd is a flag which changes the relaxivity for prohance and multihance.
% load data_info.mat
Scale=sin(data_info.FlipAngle/180*pi)/sin(data_info.ProtonDensityFlipAngle/180*pi);
% bld_T1_0=1600;  %%%% From literature
% Mz0=curves(:,1:pd_end);
% M=zeros(1,7);
T1=zeros(7,80);

 for i=1:size(curves,1)
M(i)=mean(Mz0(i,2:end-2));
 end
 M=M';
 NumberOfBaselineFrames=8;
 MaxT1=10e3;
  Curves_tmp=curves(:,pd_end+1:end);
%   TotalFrames=size(Curves_tmp,2);
  TotalFrames=size(Curves_tmp,2);
NumberOfBaselineFrames=min(NumberOfBaselineFrames,TotalFrames);
FlipAngleRatio=sin(data_info.FlipAngle/180*pi)/sin(data_info.ProtonDensityFlipAngle/180*pi);
if seriesNum==100
   SRT=data_info.SaturationRecoveryTimeAIF;
   TR=data_info.RepetitionTimeAIF;
%    Ratio=Curves_tmp./repmat(M.*FlipAngleRatio,[1 TotalFrames]);
%    Ratio=max(Ratio,1-exp(-SRT/MaxT1));
%    T1=-SRT./log(1-Ratio);
%     T1_0=-SRT./log(1-Curves_tmp./(M));
%    bld_T1_0=mean(T1(1:NumberOfBaselineFrames));
%     bld_T1_0=2000;
   n=data_info.NumberOfLinesAIF;
      
else
    SRT=data_info.SaturationRecoveryTime3D;
    TR=data_info.RepetitionTime3D;
     bld_T1_0=1530;
     n=data_info.NumberOfLines3D;
      
end
for i=1:size(Curves_tmp,2)

T1_NormalAIF=-SRT./log(1-Curves_tmp(:,i)./(M.*Scale));

for IterIndex=1:100
    Mz1=ones(size(T1_NormalAIF));
    Mxy1=zeros(size(T1_NormalAIF));
    Mz2=zeros(size(T1_NormalAIF));
    Mxy2=zeros(size(T1_NormalAIF));
     Mz2=1+(Mz2-1).*exp(-SRT./T1_NormalAIF);
%      Mz2=1+(Mz2-1).*exp(-SRT./T1);


for PulseIndex=1:n       
        Mxy1=Mxy1+Mz1.*sin(data_info.ProtonDensityFlipAngle/180*pi);
        Mz1=Mz1.*cos(data_info.ProtonDensityFlipAngle/180*pi);
        Mz1=1+(Mz1-1).*exp(-TR./T1_NormalAIF);   
%         Mz1=1+(Mz1-1).*exp(-TR./T1);
        Mxy2=Mxy2+Mz2.*sin(data_info.FlipAngle/180*pi);
        Mz2=Mz2.*cos(data_info.FlipAngle/180*pi);
         Mz2=1+(Mz2-1).*exp(-TR./T1_NormalAIF);
%         Mz2=1+(Mz2-1).*exp(-TR./T1);
end
    Mxy1=Mxy1/data_info.NumberOfLinesAIF;
    Mxy2=Mxy2/data_info.NumberOfLinesAIF;
    Scale1=sin(data_info.ProtonDensityFlipAngle/180*pi)./Mxy1;
    Mz2=1+(0-1).*exp(-SRT./T1_NormalAIF);
%     Mz2=1+(0-1).*exp(-SRT./T1);
    Scale2=Mz2.*sin(data_info.FlipAngle/180*pi)./Mxy2;
    Scale=sin(data_info.FlipAngle/180*pi)/sin(data_info.ProtonDensityFlipAngle/180*pi).*Scale1./Scale2;
    T1_NormalAIF1=-SRT./log(1-Curves_tmp(:,i)./(M.*Scale));

end
T1(:,i)=T1_NormalAIF1;
%  T1_0=repmat(mean(T1(:,1:NumberOfBaselineFrames),2),[1 TotalFrames]);
%   T1_0(i)=mean(T1(1:NumberOfBaselineFrames),2);
%%% Prohance
if flag_pd==1
r=3.8; 
else  %%%% Multihance
r=3.8;
end

%  gdcurves(:,i)=(1000./T_1(:,i)-1000./T1_0(:,i))./r;%%%% T1 is in ms
end
T1_0=repmat(mean(T1(:,1:NumberOfBaselineFrames),2),[1 TotalFrames]);
gdcurves=(1000./T1-1000./T1_0)./r;
gdcurves=abs(gdcurves);


end

  