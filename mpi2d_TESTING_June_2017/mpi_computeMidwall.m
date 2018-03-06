function mpi_computeMidwall(img,sliceNum, studyNum, outpath,referenceFrame)

[nrows, ncols, nFrames]=size(img);

imagesc(img(:,:,referenceFrame))
axis('square')
hold on
tmpImg = img(:,:,referenceFrame);

[X,Y,xCenter, yCenter, start_angle,bw1, bw2, bw3, bw5]= mpi_loadContours(img,sliceNum, studyNum, outpath,0);
myo=double(bw2)-double(bw1);
bim1=bwmorph(myo,'shrink',Inf);
figure(3)
imagesc(bim1)

[y3poly x3poly] = find(bim1);

mc = 1;
count(mc) = 1;
for i = 1:length(x3poly)-1
    if x3poly(i) == x3poly(i+1)
        count(mc) = count(mc) + 1;
        value(mc) = x3poly(i);
    else
        mc = mc + 1;
        count(mc) = 1;
        value(mc) = x3poly(i);
    end;
end;

% tmpImg2 = zeros(size(tmpImg));
% l = 0;
% for j = 1: length(count)
%     tmpImg2(y3poly(l+1):1:y3poly(l+count(j)),value(j)) = 1;
%     l = l +  count(j) ;
% end;

% -- Get the four parts of the midwall contour -- %
indexlast = [];
indexfirst = [];
k = find(count>2);
for i = 1:length(k)
    if k(i) > mean(k)
        indexlast = [indexlast k(i)];
    else
        indexfirst = [indexfirst k(i)];
    end;
end;
clear i

firstcountsum = sum(count(1:indexfirst(end)));
lastcountsum = sum(count(indexlast(1):end));

x3first(1:firstcountsum) = x3poly(1:firstcountsum);
y3first(1:firstcountsum) = y3poly(1:firstcountsum);
x3last(1:lastcountsum) = x3poly(end-(lastcountsum-1):end);
y3last(1:lastcountsum) = y3poly(end-(lastcountsum-1):end);

if mod(firstcountsum,2 ) == 0
    x3odd= x3poly(firstcountsum+1:2:end-(lastcountsum));
    y3odd = y3poly(firstcountsum+1:2:end-(lastcountsum));
    x3even = x3poly(firstcountsum+2:2:end-(lastcountsum-1));
    y3even = y3poly(firstcountsum+2:2:end-(lastcountsum-1));
else 
    x3even = x3poly(firstcountsum+1:2:end-(lastcountsum));
    y3even = y3poly(firstcountsum+1:2:end-(lastcountsum));
    x3odd = x3poly(firstcountsum+2:2:end-(lastcountsum-1));
    y3odd = y3poly(firstcountsum+2:2:end-(lastcountsum-1));
end;

[m n] = sort(y3first);
[m1 n1] = sort(y3last);

newx3poly = [];
newy3poly = [];

newx3poly = [newx3poly x3first(n(end:-1:1))];
newy3poly = [newy3poly, y3first(n(end:-1:1))];

if mod(firstcountsum,2) == 0
    newx3poly(length(newx3poly)+1:length(newx3poly)+length(x3odd)) = x3odd;
    newy3poly(length(newy3poly)+1:length(newy3poly)+length(x3odd)) = y3odd;
else 
    newx3poly(length(newx3poly)+1:length(newx3poly)+length(x3even)) = x3even;
    newy3poly(length(newy3poly)+1:length(newy3poly)+length(y3even)) = y3even;
end;

newx3poly(length(newx3poly)+1:length(newx3poly)+length(x3last)) = x3last(n1(1:end));
newy3poly(length(newy3poly)+1:length(newy3poly)+length(y3last)) = y3last(n1(1:end));

if mod(firstcountsum,2) == 0
    newx3poly(length(newx3poly)+1:length(newx3poly)+length(x3even)) = x3even(end:-1:1);
    newy3poly(length(newy3poly)+1:length(newy3poly)+length(y3even)) = y3even(end:-1:1);
else  
    newx3poly(length(newx3poly)+1:length(newx3poly)+length(x3odd)) = x3odd(end:-1:1);
    newy3poly(length(newy3poly)+1:length(newy3poly)+length(x3odd)) = y3odd(end:-1:1);
end;





newx3poly(end) = newx3poly(1);
newy3poly(end) = newy3poly(1);

bw4 = mpi_roipoly(tmpImg,newx3poly,newy3poly);
h4 = line(newx3poly, newy3poly);
set (h4, 'Color', [0.5 0.3 0.9]);

polyCoords=[newx3poly' newy3poly'];
%keyboard
tmpname=strcat(outpath,'midwall_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
save(tmpname, 'polyCoords','-ascii');

subendo = double(bw4)-double(bw1);
epi = double(bw2)-double(bw4);
figure(4)
imagesc(subendo);
colorbar

figure(5)
imagesc(epi);
colorbar

figure(6)
imagesc(subendo+epi);
colorbar

return;
