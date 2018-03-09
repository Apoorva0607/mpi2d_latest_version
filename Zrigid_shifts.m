infile =''; %% the Data structure file which has unregistered volume
load(infile);
imagesc(Unregistered(:,:,4,32));
mask=roipoly;
[RegisteredData, Offsets_XY]=jRegister3D_modified(Unregistered,mask);
[RegisteredData_new, Offsets]=jRegister3D(Unregistered,mask);
%For checking rigid shifts first run freqseolution to get a better idea of 
%where the z shift is observed, run freqseolution. Then check for manual
%shift through mpi2d stg=3.11 for that save one of the long axis file from
%freseolution as a cinemri1 file for some par file and run stg 3.11, check
%through all frames and move along slab direction to correct and the shifts
%would be stored in shiftsmanstudyxxxxslicex.txt , now next task is to read
%it.

fileID = fopen('shiftsMAN.study80000.slice8.txt','r');
data=fscanf(fileID,'%f');
data=reshape(data,[200,2]); %%% you will get shifts for all time frames along row 2 mostly.
manualcorrectionZ=0.5.*data1(2,:); %%% as the frequency has been increased twice in resolution, every shift in there corresponds to
% half the shift in time domain.
rigZ=Offsets.Z(1,:);
rigZ(rigZ>=1 & rigZ<1.5)=1;

rigZ(rigZ>=-0.5 & rigZ<0)=0;
rigZ(rigZ<0.5 & rigZ>=0)=0;
rigZ(rigZ>=0.5 & rigZ<1)=0.5;
plot(rigZ,'O')
hold on
plot(manualcorrectionZ,'+')


%%
%for line profiles
reg=squeeze(S.Registered(:,:,4,11:200));
line3=squeeze(reg(:,80,:));
figure(8);imagesc(line3);colormap gray;