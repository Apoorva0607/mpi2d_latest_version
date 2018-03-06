function out=ANTS_new_test(A,B)
i=1;
rangex=144;
rangey=144;
sx=144;
sy=144;

    filename_ref_raw=strcat('ref_',int2str(i),'.raw');
    fid = fopen (filename_ref_raw, 'wb', 'ieee-le');
    fwrite (fid, A(rangex,rangey), 'float');
    fclose(fid);
    
    filename_ref=strcat('ref_',int2str(i),'.mhd');
    fid = fopen (filename_ref, 'wb');
    header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_ref_raw);
    fprintf(fid,'%s',header_info);
    fclose(fid);
    
    
    
    filename_tar_raw=strcat('tar_',int2str(i),'.raw');
    fid = fopen (filename_tar_raw, 'wb', 'ieee-le');
    fwrite (fid, B(rangex,rangey), 'float');
    fclose(fid);
    
    filename_tar=strcat('tar_',int2str(i),'.mhd');
    fid = fopen (filename_tar, 'wb');
    header_info=sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_tar_raw);
    fprintf(fid,'%s',header_info);
    fclose(fid); 
    
    filename_out=strcat('out_',int2str(i));
    
%     filename_ref='A.tif';
%     filename_tar='B.tif';
    
    out=DR_core_new(filename_ref,filename_tar,filename_out,sx,sy,5);