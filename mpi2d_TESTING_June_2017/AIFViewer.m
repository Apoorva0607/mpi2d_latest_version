function AIFViewer()
%% things to change
slicesToAllow = [1  3 4 5 6 7 ];% 3 4 5 6 7 8];

%% adjust ratios here
% run the program and adjust the ratios here by injection amount
injectioratio = [5 5 1 2/3 2/3 2/3];
injectioratio = [1 1 1 2/3 2/3 2/3];
injectioratio = [3.6/0.7 3.6/0.7 5 1 2/3 2/3 2/3];   % for P052010A
injectioratio = [1 1 1 1 1 1 1 1];   % for P052010B

injectionsToShow = '[1 2]';
injectionsToShow = '1:length(injectionTypes)';

lineTypes = {'--*','--o','--s','--d','--^','--v','-->','--<'};


%% gather the curves based on the par files avaliable
mypwd = pwd;
if(isempty(strfind(mypwd,'Processing')))
    if(exist('Processing') == 7)
        cd('Processing');
    else
        disp('I cannot find the Processing folder.  Please naviagte to it');
        return;
    end
end
if(~isempty(strfind(mypwd,'Output')))
    cd('..');
end
parfiles = dir('*.par');
curveFileNames = {};
clear temp; temp.hello = 0;
mymap = containers.Map({5},{temp});
remove(mymap,5);

for p=1:length(parfiles)
    if(strcmp(parfiles(p).name,'mpi2d.par')) continue; end    %kill off mpi2d
    if(~isempty(strfind(parfiles(p).name,'rad'))) continue; end   %kill off subsets
    
    %extract the series and slice number
    A = sscanf(parfiles(p).name,'series%d.slice%d.par');
    if(isempty(A)) continue; end
    series = A(1); slice = A(2); %subset = A(3);
    if(isempty(find(slicesToAllow == slice)))
        continue;
    end
    
    %get at the injection type
    ParFileName = parfiles(p).name;
    ReadPar
    where = strfind(infile,'/');
    dicomFolder = infile((where(end-1)+1):(where(end)-1));
    parenStart = strfind(dicomFolder,'(');
    injectionType = dicomFolder(1:(parenStart-1));
    
    %find the delteSIcurves file
    template = strrep(strrep(parfiles(p).name,'series','study'),'par','mat');
    %curvesCandidates = dir(['Output/gdcurves*' template]);
    curvesCandidates = ['gdcurves.study',int2str(series),'.slice',int2str(slice),'.AIF_',int2str(seriesNumAIF),'_',int2str(1),'_1.mat'];
    %throw it in the map
    if(~isKey(mymap,injectionType))   %DEV
        clear temp
        load(['Output/' curvesCandidates])
        temp.curves = containers.Map({[0]},{gdcurves});
        temp.series = containers.Map({[0]},{[0]});
        remove(temp.curves,[0]);
        mymap(series) = temp;
    end
    data = mymap(series);
    data.curvesFileName = curvesCandidates;
    load(['Output/' curvesCandidates]);
    data.curves(slice) = gdcurves;
    data.series(slice) = series;
    
end

disp('----------------slicei--');
%% display the curves
figure(1);
%set(gcf,'Color',[.23,0.44,0.34])
hold on
legendStrings = {};
%mycolormap = lines;    % change this to be whatever colormap you want the curves to divide between themselves
mycolormap = hsv;   %the curves are evenly spaced using the colormap given here
%mycolormap = gray;   %use this one if you want to print things out and kill the 'Color',[.23 ,.44...] thing above 
injectionTypes = keys(mymap);
for injectioni = eval(injectionsToShow)
    injection = injectionTypes{injectioni};
    disp(['Ratio(' injection ') = ' num2str(injectioratio(injectioni))])
    
    data = mymap(injection);
    slices = keys(data.curves);
    for slicei = 1:length(slices)
        slice = slices{slicei};
        curves = data.curves(slice);
        AIF = curves(8,:);
                
        offset=0;
        AIF(1:max(1,offset)) = 0;
        AIF = circshift(AIF,[0 offset]);
        AIF(1:max(1,offset)) = 0;
        AIF = AIF*injectioratio(injectioni);
        plot(AIF,lineTypes{injectioni},'Color',1-mycolormap(1+floor((injectioni-1)*length(mycolormap)/length(injectionTypes)),:),'LineWidth',.5);
        legendStrings{length(legendStrings)+1} = ['(' num2str(data.series(slice)) ',' num2str(slice) ')'];
    end
end
hold off
xlabel('Time Frame');
ylabel('[Gd]');
title({'AIF seperated by color and marker type'})
hold off
if(isempty(legendStrings)) 
    return; 
end
lineToExpress = ['legend(''' strrep(legendStrings{1},'_','\_') ''''];
for line=2:length(legendStrings)
    lineToExpress = [lineToExpress ',' '''' strrep(legendStrings{line},'_','\_') ''''];
end
lineToExpress = [lineToExpress ');'];
eval(lineToExpress);
end