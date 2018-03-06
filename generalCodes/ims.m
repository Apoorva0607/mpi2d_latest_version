function h = ims(varargin)

f = [];
img = [];
indexedImage = false;
rgbImage = false;
rng = [];

[img,args,nvpargs] = parseArguments(varargin);

for i=1:length(args)
    if (isa(args{i},'function_handle'))
        f = args{i};
        args = {args{1:i-1} args{i+1:length(args)}};
        break;
    end;
end;

% identity operator
if (isempty(f))
    f = @(x)(x);
end;

img = squeeze(f(img));

if (ndims(img) == 2)
    indexedImage = true;
elseif (ndims(img) == 3 && size(img,3) == 3)
    rgbImage = true;
end;

if (~indexedImage && ~rgbImage); 
    error('ims :: invalid image data');
end;

args = {img args{:}};
h = imagesc(args{:});

if (isfield(nvpargs,'Box'))
    rng = getfield(nvpargs,'Box');
   
    if (~iscell(rng)) rng = {rng}; end;
end;

% draw box
for i=1:length(rng)
    if (~isempty(rng{i}))
        sz = size(img);
                    
        xlo = rng{i}(2,1);
        xhi = rng{i}(2,2);
        ylo = rng{i}(1,1);
        yhi = rng{i}(1,2);

        hold on;
        line([xlo xhi],[ylo ylo],'Color','w');
        line([xlo xhi],[yhi yhi],'Color','w');
        line([xlo xlo],[ylo yhi],'Color','w');
        line([xhi xhi],[ylo yhi],'Color','w');
        hold off;
    end;
end;

return;
