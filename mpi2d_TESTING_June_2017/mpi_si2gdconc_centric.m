function gd_curves=mpi_si2gdconc_centric(curves,sliceNum,studyNum, outpath, numPreContrast_bld, numPreContrast_tiss, numSkip, scale_to_study1,ranget)

% changing this 11/13/06 to more analytical. Calling old version
% mpi_si2gdconcByLUT  (that version was based on imaging vials of known
% concentrations

if nargin==0,
    error('sliceNum argument needed for mpi_si2gdconc.');
end




figure(2)
showcurves([bldGds(:)'; tissGds(:,:)'],'[Gd] (mmol)');
title('From analytical centric expression')
%pause

gd_curves =[ bldGds'; tissGds'];  % so in same format as curves

return;
