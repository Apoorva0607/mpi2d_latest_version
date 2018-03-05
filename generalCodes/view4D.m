% Interactively explore 4D data sets 
%   Displays image of current spatial slice that can be changed by scroll wheel
%   and time series corresponding to current mouse position. If multiple input
%   4D data sets are provided (via cell array of 4D arrays), the first array is
%   used for image data and time curves are plotted for all data sets.
%
%   Example usage : view4D({data1,data2,data3},imageRange)
%
%   Options :
%       'LineColormap'  - use specified colormap for plotting lines
%       'ImageColormap' - use specified colormap for displaying images
% 
%   Copyright 2011 Matthias Christian Schabel (matthias @ stanfordalumni . org)
%   Oregon Health & Science University
%   Advanced Imaging Research Center
%   3181 SW Sam Jackson Park Road L452
%   Portland, OR 97239
%   
function view4D(varargin)

[alldata,args,nvpargs] = parseArguments(varargin);

if (~iscell(alldata)) alldata = {alldata}; end;

numberOfDataSets = length(alldata);

try
    lineColormap = nvpargs.LineColormap;
catch e
    lineColormap = prism(numberOfDataSets);
    lineColormap(1,:) = [0 0 0];
end;

try
    imageColormap = nvpargs.ImageColormap;
catch e
    imageColormap = colormap;
end;

data = alldata{1};

for i=1:numberOfDataSets
    if (ndims(alldata{i}) ~= 4) error('DCELAB::fourDExplorer : data must be 4D'); end;
    if (~isequal(size(alldata{i}),size(data))) error('DCELAB::fourDExplorer : data sets are different sizes'); end;
end;

sz = size(data);

timePoints = sz(1);
numberOfSlices = sz(2);
imsz = sz(3:4);

timeIdx = ceil(timePoints/10);
lasttimeIdx = 0;

mn = nanmean(reshape(data,[sz(1) prod(sz(2:end))]),2);

%keyboard
%timeIdx = first(find(mn == nanmax(mn)));
timeIdx = find(mn == nanmax(mn));

zidx = ceil(numberOfSlices/2);
xidx = 1;
yidx = 1;
lastzidx = 0;
lastxidx = 0; 
lastyidx = 0;

f1 = figure;
set(f1,'Visible','off');
set(f1,'Position',[50 50 600 840]);

h = [];
ln = [];

ax1 = axes('Position',[.05 .05 .9 .45]);
plotImage(ax1,timeIdx,zidx);
colormap(ax1,imageColormap);

ax2 = axes('Position',[.05 .55 .9 .4]);
plotTimeSeries(ax2,zidx,xidx,yidx);
set(f1,'Visible','on');

% set(f1,'WindowScrollWheelFcn',@scrollFigure);
set(f1,'WindowButtonMotionFcn',@mouseMoved);

    function mouseMoved(src,evnt)            
        currPoint = get(ax1,'CurrentPoint');

        xidx = floor(currPoint(1,2)+.5);
        yidx = floor(currPoint(1,1)+.5);

        % can only scroll if more than one slice
        if (numberOfSlices > 1)
            set(f1,'WindowScrollWheelFcn',@scrollFigure);
        else
            set(f1,'WindowScrollWheelFcn',[]);
        end;

        % can only scroll in time if more than one time point
        if (timePoints > 1)
            set(f1,'WindowButtonDownFcn',@startDrag);
            set(f1,'WindowButtonUpFcn',@stopDrag);
        else
            set(f1,'WindowButtonDownFcn',[]);
            set(f1,'WindowButtonUpFcn',[]);
        end;
        
        plotTimeSeries(ax2,zidx,xidx,yidx);
    end

    function scrollFigure(src,evnt)
        if evnt.VerticalScrollCount > 0
            if (zidx < numberOfSlices)
                zidx = zidx+1;
            end;
        elseif evnt.VerticalScrollCount < 0 
            if (zidx > 1)
                zidx = zidx-1;
            end;
        end

        lastzidx = 0;
        plotImage(ax1,timeIdx,zidx);
        plotTimeSeries(ax2,zidx,xidx,yidx);
    end 

    function plotImage(ax,timeIdx,zidx)
        axes(ax);
        newargs = {data(timeIdx,zidx,:,:) args{:}};
        h = ims(newargs{:});
        axis equal;
        set(gca,'XTick',[],'YTick',[])        
        title(['Slice ' num2str(zidx)]);
        keyboard
        drawnow;
    end 

    function plotTimeSeries(ax,zidx,xidx,yidx)        
        % only plot if x and y are within image and pixel position has changed
        if (xidx >= 1 && xidx <= imsz(1) && yidx >= 1 && yidx <= imsz(2) && (xidx ~= lastxidx || yidx ~= lastyidx || zidx ~= lastzidx || timeIdx ~= lasttimeIdx))
            lasttimeIdx = timeIdx;
            
            lastzidx = zidx;
            lastxidx = xidx;
            lastyidx = yidx;

            curves = zeros([numberOfDataSets sz(1)]);           
             
            for i=1:numberOfDataSets
                curves(i,:) = squeeze(alldata{i}(:,zidx,xidx,yidx));
            end;
           
            h = plot(ax,1:timePoints,curves,...
                        timeIdx,squeeze(data(timeIdx,zidx,xidx,yidx)),'ro');
            
            set(h(1),'LineStyle','-','LineWidth',1,'Marker','.','Color',[0 0 0]);
            
            for i=2:numberOfDataSets
                set(h(i),'LineStyle','-','LineWidth',.5,'Color',lineColormap(i,:));
            end;
            
            xlim = get(ax,'XLim');
            ylim = get(ax,'YLim');
            
            ln = line('Parent',ax,'XData',[timeIdx timeIdx],'YData',ylim,'Color','r');
            
            val = data(timeIdx,zidx,xidx,yidx);
            text('Parent',ax,'Position',[9 .5],'Units','centimeters',...
                 'String',['(' num2str(timeIdx) ' | ' num2str(zidx) ',' num2str(xidx) ',' num2str(yidx) ')'  ' : ' num2str(val)]);
            drawnow;
        end
    end

    function startDrag(src,evnt)         
        set(f1,'WindowButtonMotionFcn',{@mouseDrag,get(ax2,'CurrentPoint'),timeIdx});
        set(f1,'WindowButtonUpFcn',@stopDrag);
    end

    function mouseDrag(src,evnt,startPt,startTimeIdx)       
        currPt = get(ax2,'CurrentPoint');
       
        timeIdx = bound(round(currPt(1,1)-startPt(1,1))+startTimeIdx,1,timePoints);
        
        set(ln,'XData',[timeIdx timeIdx]);
        
        plotImage(ax1,timeIdx,zidx);
        plotTimeSeries(ax2,zidx,xidx,yidx);
    end

    function stopDrag(src,evnt)        
        currPoint = get(ax1,'CurrentPoint');

        xidx = floor(currPoint(1,2)+.5);
        yidx = floor(currPoint(1,1)+.5);
        
        plotTimeSeries(ax2,zidx,xidx,yidx);

        set(f1,'WindowButtonMotionFcn',@mouseMoved);
    end

end % scroll_wheel
