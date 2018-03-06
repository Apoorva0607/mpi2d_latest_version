function copy_process(parent_ser,parent_slice,seriesNum,sliceNum)

load(['Output/blood_polyCoords.study',int2str(parent_ser),'.slice',int2str(parent_slice)]);

load(['Output/endo_polyCoords.study',int2str(parent_ser),'.slice',int2str(parent_slice)]);

load(['Output/epi_polyCoords.study',int2str(parent_ser),'.slice',int2str(parent_slice)]);

load(['Output/Roi_start_angle.study',int2str(parent_ser),'.slice',int2str(parent_slice)]);

load(['Output/shifts.study',int2str(parent_ser),'.slice',int2str(parent_slice),'.txt']);


f_tar=['Output/blood_polyCoords.study',int2str(seriesNum),'.slice',int2str(sliceNum)];
tmp=['blood_polyCoords_study',int2str(seriesNum),'=','blood_polyCoords_study',int2str(parent_ser)];
var_tar=[' blood_polyCoords_study',int2str(seriesNum)];
eval(tmp);
tmp1=['save ',f_tar,var_tar,' -ascii'];
eval(tmp1)

f_tar=['Output/endo_polyCoords.study',int2str(seriesNum),'.slice',int2str(sliceNum)];
tmp=['endo_polyCoords_study',int2str(seriesNum),'=','endo_polyCoords_study',int2str(parent_ser)];
var_tar=[' endo_polyCoords_study',int2str(seriesNum)];
eval(tmp);
tmp1=['save ',f_tar,var_tar,' -ascii'];
eval(tmp1)

f_tar=['Output/epi_polyCoords.study',int2str(seriesNum),'.slice',int2str(sliceNum)];
tmp=['epi_polyCoords_study',int2str(seriesNum),'=','epi_polyCoords_study',int2str(parent_ser)];
var_tar=[' epi_polyCoords_study',int2str(seriesNum)];
eval(tmp);
tmp1=['save ',f_tar,var_tar,' -ascii'];
eval(tmp1)

f_tar=['Output/Roi_start_angle.study',int2str(seriesNum),'.slice',int2str(sliceNum)];
tmp=['Roi_start_angle_study',int2str(seriesNum),'=','Roi_start_angle_study',int2str(parent_ser)];
var_tar=[' Roi_start_angle_study',int2str(seriesNum)];
eval(tmp);
tmp1=['save ',f_tar,var_tar,' -ascii'];
eval(tmp1)

f_tar=['Output/shifts.study',int2str(seriesNum),'.slice',int2str(sliceNum),'.txt'];
tmp=['shifts_study',int2str(seriesNum),'_slice',int2str(sliceNum),'=','shifts_study',int2str(parent_ser),'_slice',int2str(parent_slice)];
var_tar=[' shifts_study',int2str(seriesNum),'_slice',int2str(sliceNum),];
eval(tmp);
tmp1=['save ',f_tar,var_tar,' -ascii'];
eval(tmp1)