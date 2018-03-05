load TestData;
data_info.ContrastType=1;
data_info.NumberOfBaselineFrames=1;

% Find T1 using my method
% gdcurves is not correct since this is not a dynamic data set
[gdcurves T1]=conc_gd_v2(curves,100,data_info);

% Sort vials from smallest to largest [Gd]
T1=T1([1 9 2 3 4 8 7 6 5 10 11 12 13]);

% First two vials are the saline (no Gd)
T10=mean(T1(1:2));

% Get [Gd]
r=3.8;
C1=-(1000./T10-1000./T1)./r;

% Find T1 using equation
% gdcurves is not correct since this is not a dynamic data set
[gdcurves1 T1]=conc_gd_v3(curves,100,data_info);

% Sort vials from smallest to largest [Gd]
T1=T1([1 9 2 3 4 8 7 6 5 10 11 12 13]);

% First two vials are the saline (no Gd)
T10=mean(T1(1:2));

% Get [Gd]
r=3.8;
C2=-(1000./T10-1000./T1)./r;