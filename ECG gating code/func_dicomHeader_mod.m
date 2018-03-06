function func_dicomHeader_mod(a,templatedirtxt, indata, outname, seriesnum, seriesdesc); 

addDicomHeader(templateFilename, indata(:,:,k), seriesnum,seriesdesc, k, outname);