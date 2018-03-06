% ----------------------------------------------------------------------
% [gdcurves T1]=conc_gd_v3(curves,seriesNum,data_info)
% UCAIR - June 2017
% ----------------------------------------------------------------------
% Converts signal into Gd concentration.
%
% The input argument curves should be size [NumberOfRegions,TotalFrames]
% where the first frames are proton density data. The number of proton
% density frames is determined by data_info. It is also necessary to
% specify the number of baseline frames (time points before contrast
% injection).
%
% The input argument seriesNum determines if this is 2D AIF data or 3D
% tissue data. Currently seriesNum 100 and 101 correspond to the 2D AIF.
%
% The input 'data_info' should contain the following fields:
% ContrastType: 1 for ProHance (default), 0 for MultiHance
% NumberOfBaselineFrames: Frames before injection (excluding proton density)
% NumberOfProtonDensityFrames: Proton density measurements
% FlipAngle: Flip angle of tissue data (degrees)
% ProtonDensityFlipAngle: Flip angle of proton density data (degrees)
% SaturationRecoveryTimeAIF: For AIF only (ms)
% RepetitionTimeAIF: For AIF only (ms)
% NumberOfLinesAIF: For AIF only
% SaturationRecoveryTime3D: For 3D only (ms)
% RepetitionTime3D: For 3D only (ms)
% LinesInCentralPartition: For 3D only (ms)
%
% If you assume perfect saturation recovery, then the T1 relaxation time
% can be calculated from the ratio of proton density and tissue
% longitudinal magnetizations: MzTissue/MzProtonDensity=1-exp(-SRT/T1).
% However, we only measure the transverse magnetization. If the proton
% density and tissue data are acquired with the same excitation flip angle
% then the ratio of the transverse magnetizations is the same as the ratio
% of the longitudinal magnetizations and no correction is necessary.
% However, if the proton density and the tissue data are acquired with
% different flip angles then each transverse signal must be corrected to
% account for the varying flip angles. This is done by dividing each
% transverse signal by sin() of the corresponding flip angle. Since we are
% only interested in the ratio of the two transverse signals, the two
% correction factors can be combined into a single correction factor which
% we call 'FlipAngleRatio'.
%
% This correction factor assumes there is a single RF exciation after each
% saturation recovery period. If there are multiple RF pulses, then these
% must also be taken into consideration. In this case we calculate
% additional correction factors to account for the difference between a
% single and multiple RF pulses.
% ----------------------------------------------------------------------
function [gdcurves T1]=conc_gd_v3(curves,seriesNum,data_info)

% Setup for conversion
if (isfield(data_info,'ContrastType'))
    ContrastType=data_info.ContrastType;
else
    ContrastType=1;
end
if (isfield(data_info,'NumberOfBaselineFrames'))
    NumberOfBaselineFrames=data_info.NumberOfBaselineFrames;
else
    NumberOfBaselineFrames=1;
end
if (isfield(data_info,'NumberOfProtonDensityScans'))
    NumberOfProtonDensityFrames=data_info.NumberOfProtonDensityScans;
else
    NumberOfProtonDensityFrames=1;
end
FlipAngleRatio=sin(data_info.FlipAngle/180*pi)/sin(data_info.ProtonDensityFlipAngle/180*pi);
NumberOfIterations=100;
MaxT1=10e3; % This is used to prevent divide by zero in the log function
if ((seriesNum==100)||(seriesNum==101)) 
    SRT=data_info.SaturationRecoveryTimeAIF;
    TR=data_info.RepetitionTimeAIF;
    NumberOfPulses=data_info.NumberOfLinesAIF;
else
    SRT=data_info.SaturationRecoveryTime3D;
    TR=data_info.RepetitionTime3D;
    NumberOfPulses=data_info.LinesInCentralPartition;
end

% Define signals
% The first time frames are the proton densty signal
% Ignore the first proton density time frame as we are not in steady state (if possible)
% Ignore the last proton density time frame due to temporal blurring with TCR (if possible)
% Ignore the first tissue time frame due to emporal blurring with TCR (if possible)
if (NumberOfProtonDensityFrames>2)
    ProtonDensitySignal=mean(curves(:,2:(NumberOfProtonDensityFrames-1)),2);
else
    ProtonDensitySignal=mean(curves(:,1:NumberOfProtonDensityFrames),2);
end
TissueSignal=curves(:,min((NumberOfProtonDensityFrames+2),size(curves,2)):size(curves,2));
TotalFrames=size(TissueSignal,2);
NumberOfBaselineFrames=min(NumberOfBaselineFrames,TotalFrames);
    
% First T1 estimate assumes a single RF pulses
% MzTissue/MzProtonDensity=1-exp(-SRT/T1)
% MzTissue/MzProtonDensity=MxyTissue/(MxyProtonDensity*FlipAngleRatio) 
% MxyTissue/(MxyProtonDensity*FlipAngleRatio)=1-exp(-SRT/T1)
% See function comments above for more details
Ratio=TissueSignal./repmat(ProtonDensitySignal.*FlipAngleRatio,[1 TotalFrames]);
Ratio=max(Ratio,1-exp(-SRT/MaxT1));
T1=-SRT./log(1-Ratio);

% Estimate T1 for multiple RF pulses with an iterative method
for IterIndex=1:NumberOfIterations
    % Simulate signal evolution with current T1 estimate
    % Normalize the equilibrium magnetization (M0=1.0)
    MzTissue=ones(size(TissueSignal));
    MxyTissue=zeros(size(TissueSignal));
    MzProtonDensity=ones(size(ProtonDensitySignal));
    MxyProtonDensity=zeros(size(ProtonDensitySignal));
    PreContrastT1=mean(T1(:,1:NumberOfBaselineFrames),2); % T1 before injection
    % Apply saturation pulse to tissue signal and recover during SRT
    % No saturation pulse is applied during proton density scans
    % Simulate the effect of each RF pulse and the recovery during each TR
    % Proton density recovers with PreContrastT1
    % Assume image contrast is an average of transverse signal over all pulses
    % Use equation from Jerosch-Herold et al., Med Phys 1998;25:73-84.
    for PulseIndex=1:NumberOfPulses
        Alpha=cos(data_info.FlipAngle/180*pi)*exp(-TR./T1);
        MzTissue=(1-exp(-SRT./T1)).*(Alpha.^(PulseIndex-1))+(1-exp(-TR./T1)).*(1-Alpha.^(PulseIndex-1))./(1-Alpha); % Longitudinal change due to RF
        MxyTissue=MxyTissue+MzTissue.*sin(data_info.FlipAngle/180*pi); % Excited transverse signal
        Alpha=cos(data_info.ProtonDensityFlipAngle/180*pi)*exp(-TR./T1);
        MzProtonDensity=1.*(Alpha.^(PulseIndex-1))+(1-exp(-TR./T1)).*(1-Alpha.^(PulseIndex-1))./(1-Alpha); % Longitudinal change due to RF
        MxyProtonDensity=MxyProtonDensity+MzProtonDensity.*sin(data_info.ProtonDensityFlipAngle/180*pi); % Excited transverse signal
    end
    MxyTissue=MxyTissue./NumberOfPulses; % Kspace center is sampled each RF pulse so average
    MxyProtonDensity=MxyProtonDensity./NumberOfPulses; % Kspace center is sampled each RF pulse so average
    % Correction for Proton Density Signal
    % MxyProtonDensity is the predicted transverse signal accounting for multiple RF pulses after the SRT
    % Mxy after a single RF pulse folowing the SRT is: M0*sin(FlipAngle)
    % M0=1.0, no saturation, then Mxy=sin(data_info.ProtonDensityFlipAngle/180*pi)
    % The correction between the single RF and the multiple RF is:
    % ProtonDensityScale=Mxy/MxyProtonDensity
    ProtonDensityScale=sin(data_info.ProtonDensityFlipAngle/180*pi)./MxyProtonDensity;
    % Correction for Tissue Signal
    % MxyTissue is the predicted transverse signal accounting for multiple RF pulses after the SRT
    % Mxy after a single RF pulse is: M0*(1-exp(-SRT/T1))*sin(FlipAngle)
    % M0=1.0, no saturation, then Mxy=(1-exp(-SRT/T1)).*sin(data_info.FlipAngle/180*pi)
    % The correction between the single RF and the multiple RF is:
    % TissueScale=Mxy/MxyTissue
    TissueScale=(1-exp(-SRT./T1)).*sin(data_info.FlipAngle/180*pi)./MxyTissue;
    % Find the new estimate of T1
    Ratio=(TissueSignal.*TissueScale)./repmat((ProtonDensitySignal.*ProtonDensityScale).*FlipAngleRatio,[1 TotalFrames]);
    Ratio=max(Ratio,1-exp(-SRT/MaxT1));
    T1=-SRT./log(1-Ratio);
end

% Convert from T1 to [Gd]
if (ContrastType==1) % ProHance
    r=2.8; 
else % MultiHance
    r=3.8;
end
 T1_0=repmat(mean(T1(:,1:NumberOfBaselineFrames),2),[1 TotalFrames]);
% T1_0=1800;
gdcurves=(1000./T1-1000./T1_0)./r; % T1 is in ms 