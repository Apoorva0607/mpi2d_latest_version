function [coeff] = regressPlot_ed(X,Y,offset)


%y=mx+b
tmp=polyfit(X,Y,1);
m=tmp(1); b=tmp(2)

r=corrcoef(Y,X);
disp(['Fit eqtn: y=',num2str(m),'x+',num2str(b),' r=',num2str(r(1,2))]);

	

minx = min(X);
maxx = 1.2*max(X);
maxx=max(maxx,1.1*max(Y));
t1 = linspace(minx,maxx,16);
t1 = linspace(0,maxx,36);
XI = [ones(1,36);t1];
YI = m*t1 + b;

%plot(X+offset-1,Y,'g*','linewidth',2); hold on
plot(t1+offset-1,YI(1,:),'k','linewidth',2.5)
xlim([0 maxx]); ylim([0 maxx])
axis square
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',14)
% hmmm, line seems a little crooked at top...
% evd changed to just plot out regression of Y on X
%plot(X,Y,'+',XI(2,:),YI(1,:),XI(2,:),YI(2,:),':')

%%%% to automatically display the linear regression eqns and R^2 values on the plot 
xx=32.8;yy=3.8; %location of eqn displayed on plot
%xx=32.8;yy=4.8; %location of eqn displayed on plot

xx=1.02;yy=0.35; %location of eqn displayed on plot
labelString1=sprintf('y=%1.2f*x+%1.2f',m,b);
text(xx,yy,labelString1,'FontSize',14)
% xx=32.8;yy=3.5  %0.5;
% %xx=32.8;yy=4.5  %0.5;
xx=1.02; yy=0.15
 labelString2=sprintf('r=%1.2f',r(1,2) );
 text(xx,yy,labelString2,'FontSize',14)


return
