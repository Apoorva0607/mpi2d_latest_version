function [coeff] = regressPlot_ed_bmi(X,Y,offset, bmi)


%y=mx+b
tmp=polyfit(X,Y,1);
m=tmp(1); b=tmp(2)

r=corrcoef(Y,X);
disp(['Fit eqtn: y=',num2str(m),'x+',num2str(b),' r=',num2str(r(1,2))]);

	

minx = min(X);
maxx = 1.2*max(X);
maxx=1.1*max(Y);
t1 = linspace(minx,maxx,16);
t1 = linspace(0,maxx,16);
XI = [ones(1,16);t1];
YI = m*t1 + b;

plot(X+offset-1,Y,'g*','linewidth',2); hold on
plot(t1+offset-1,YI(1,:),'k','linewidth',2.5)
xlim([0 maxx]); ylim([0 maxx])
axis square
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
% hmmm, line seems a little crooked at top...
% evd changed to just plot out regression of Y on X
%plot(X,Y,'+',XI(2,:),YI(1,:),XI(2,:),YI(2,:),':')

%%%% to automatically display the linear regression eqns and R^2 values on the plot 
xx=2.8;yy=1.2; %location of eqn displayed on plot
labelString1=sprintf('y=%1.2f*x+%1.2f',m,b);
text(xx,yy,labelString1,'FontSize',14)
xx=2.8;yy=0.5;
labelString2=sprintf('r=%1.2f',r(1,2) );
text(xx,yy,labelString2,'FontSize',14)



figure
% try plot with color/intensity mapped to weight, 
load spect.cmap
spect=spect./max(spect(:));  % needs to be 0 to 1
% plot(x,y,'o','LineWidth',2,'Color',spect(201,:))

bmiNorm=(bmi-min(bmi))./(max(bmi)-min(bmi));  % want max to be 1 and min to be 0 to use whole range. 
bmiNorm=(bmi-10)./(max(bmi)-10);  % actually bottom of range all looks the same, so don't take out min... well, take out 5 or 10.. change set(hh) below to match



for i=1:length(X)
    plot(X(i)+offset-1,Y(i),'*','linewidth',3.9,'Color',[0,bmiNorm(i),0]); hold on
%     tmpp=round(bmiNorm(i)*240)+1
%     plot(X(i)+offset-1,Y(i),'*','linewidth',3,'Color',spect(tmpp,:)); hold on
end
plot(t1+offset-1,YI(1,:),'k','linewidth',2.5)
xlim([0 maxx]); ylim([0 maxx])
axis square
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)

cc=[0:.01:1;zeros(1,101);zeros(1,101)]';
cc2(:,2)=cc(:,1); cc2(:,1)=cc(:,2); cc2(:,3)=cc(:,3);
 
colormap(cc2)
%colormap(spect)
colorbar
% change label for max and min of color:
hh=gca;
%set(hh,'Clim',[floor(min(bmi)), max(bmi)])
set(hh,'Clim',[10, max(bmi)])
% set(hh,'Clim',[floor(min(bmi)), (255/240)*max(bmi)])
% set(hh,'Clim',[floor(min(bmi)), max(bmi)])

%%%% to automatically display the linear regression eqns and R^2 values on the plot 
xx=2.8;yy=1.2; %location of eqn displayed on plot
labelString1=sprintf('y=%1.2f*x %1.2f',m,b);
text(xx,yy,labelString1,'FontSize',14)
xx=2.8;yy=0.5;
labelString2=sprintf('r=%1.2f',r(1,2) );
text(xx,yy,labelString2,'FontSize',14)


return
