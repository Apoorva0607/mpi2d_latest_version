function T1map_dicoms(T1,inpath,instring,outpath,seriesNum,seriesDesc)


cd('Pre_PD_SA_40')

list=dir('IM*');
 cd ..

for i=1:length(list)
    
%  cd(list(i).name)
outname=['Dicom/slice',num2str(i),'.dcm'];
addDicomHeader(list(i).name,(T1(:,:,i)),seriesNum,seriesDesc,1,outname)





end

