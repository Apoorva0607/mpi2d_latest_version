

flowDir='/v/raid1a/apedgaon/MRIdata/Cardiac/Prisma_version1/P102516/Processing/Output/';
cd(flowDir)
AIF_series=100;
AIF_slice=1;
study=[3000,4000,5000,6000,7000];
slices=[3,4,5,6,7];
[ ktrans_sys_stress] = findFlowFiles(flowDir, study, slices,AIF_series,AIF_slice);

% study=[13000,14000,15000,16000,17000,18000];
% slices=[3,4,5,6,7,8];
% [ ktrans_dias_rest] = findFlowFiles(flowDir, study, slices,AIF_series,AIF_slice);

cd ..
 save('ktrans_sys_stress.mat','ktrans_sys_stress');
%  save('ktrans_dias_rest.mat','ktrans_dias_rest');

