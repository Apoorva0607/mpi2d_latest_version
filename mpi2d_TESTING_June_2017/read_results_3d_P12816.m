

flowDir='/v/raid1a/apedgaon/MRIdata/Cardiac/Prisma_version1/P111716/Processing_stress/Output/';
cd(flowDir)
AIF_series=100;
AIF_slice=1;
study=[4000,5000,6000,7000,8000];
slices=[4,5,6,7,8];
[ ktrans_sys_rest] = findFlowFiles(flowDir, study, slices,AIF_series,AIF_slice);

% study=[13000,14000,15000,16000,17000,18000];
% slices=[3,4,5,6,7,8];
% [ ktrans_dias_rest] = findFlowFiles(flowDir, study, slices,AIF_series,AIF_slice);

cd ..
 save('ktrans_sys_rest.mat','ktrans_sys_rest');
%  save('ktrans_dias_rest.mat','ktrans_dias_rest');

