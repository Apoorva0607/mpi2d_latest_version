function Vb = getVb(f)
fid=fopen(f,'r');
if(fid < 0)
    disp(['Couldn''t open file: ' f]);
    Vb = -1;
    return
end
Vb = [];
tline = fgetl(fid);
while ischar(tline)
    if(~isempty(strfind(tline,'fv  ')))
        A = sscanf(tline,'fv  %d     %f');
        Vb(A(1)) = A(2);
    end
    tline = fgetl(fid);
end

fclose(fid);