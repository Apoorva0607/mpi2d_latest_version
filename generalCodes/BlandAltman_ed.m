function [means, diffs, meanDiff, confRang] = BlandAltman_ed(dat1,dat2,weights)


fPlot=1; % default plot


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
% if weights>cutoff
%         diffs1=dat1-dat2;
%     else
%         diffs2=
% end
meanDiff = mean(diffs);
 stdDiff_old = std(diffs);

[p,tab]=anova1(reshape(diffs,[18,10]));

var1=(tab{2,4}-tab{3,4})/18;
var2=tab{3,4};
var_act=var1+var2;
stdDiff=sqrt(var_act);
confRang = [meanDiff + 1.96*stdDiff, meanDiff - 1.96*stdDiff]
confRang_old = [meanDiff + 1.96*stdDiff_old, meanDiff - 1.96*stdDiff_old]

if(fPlot)
    figure(2)
    minMeans=min(means);
    maxMeans=max(means);
    x=minMeans-0.1:0.1:maxMeans+0.1;
    clf;
    
     plot(means(1:18),diffs(1:18),'ro')
     hold on
    plot(means(19:36),diffs(19:36),'go');
    
    plot(means(37:54),diffs(37:54),'bo');
%    
   plot(means(55:72),diffs(55:72),'mo');
   plot(means(73:90),diffs(73:90),'co');
   plot(means(91:108),diffs(91:108),'ko');
   plot(means(109:126),diffs(109:126),'yo');
   plot(means(127:144),diffs(127:144),'Marker','o','Color',[1 0.5 0.4]);
   plot(means(145:162),diffs(145:162),'Marker','o','Color',[0.5 0.5 1]);
   plot(means(163:180),diffs(145:162),'Marker','o','Color',[0.7 0.3 0.8]);
    hold on
    plot(x,ones(1,length(x)).*confRang(1),'k-.','linewidth',2)
    plot(x,ones(1,length(x)).*confRang(2),'k-.','linewidth',2)
    plot(x,ones(1,length(x)).*meanDiff,'k:','linewidth',2)
    set(gcf,'Color',[1 1 1])
end