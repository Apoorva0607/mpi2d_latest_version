function [img,args,nvpargs] = parseArguments(argin)

args = {};
nvpargs = struct();

i = 1;
j = 1;

while i<=length(argin)
    if (isa(argin{i},'char'))
        nvpargs = setfield(nvpargs,argin{i},argin{i+1});
        i = i+2;
    else
        args{j} = argin{i};
        i = i+1;
        j = j+1;
    end;
end;

for i=1:length(args)
    if (ndims(squeeze(args{i})) >= 2)
        img = args{i};
        args = {args{1:i-1} args{i+1:length(args)}};
        break;
    end;  
end;

return;

