function [means, diffs, meanDiff, confRang] = BlandAltman(dat1,dat2,fPlot)

if nargin < 3
    fPlot=1; % default plot
end

if nargin <2
    error('BlandAltman :: insufficient data')
end

if size(dat1,1) ~= 1
    dat1=dat1';
end
if size(dat2,1) ~= 1
    dat2=dat2';
end

means = mean([dat1; dat2]);
diffs = dat1-dat2;

meanDiff = mean(diffs);
stdDiff = std(diffs);

confRang = [meanDiff + 1.96*stdDiff, meanDiff - 1.96*stdDiff];

if(fPlot)
    minMeans=min(means);
    maxMeans=max(means);
    x=minMeans-0.1:0.1:maxMeans+0.1;
    clf;
    plot(means,diffs,'ro')
    hold on
    plot(x,ones(1,length(x)).*confRang(1),'k-.')
    plot(x,ones(1,length(x)).*confRang(2),'k-.')
    plot(x,ones(1,length(x)).*meanDiff,'k:')
    hold off
end