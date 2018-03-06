function [coeff] = twoplot(Y,X)

% twoplot    Plot regression of Y on X and X on Y.
%
%            Usage:
%
%            [coeff] = twoplot(Y,X)
%
%            Returns:
%
%            coeff  - column 1 gives coefficients of regression Y on X
%                     column 2 gives rearranged coefficients of X on Y

% M J Chlond - May 96
% m.chlond@uclan.ac.uk

% modified a little to print out t-values of slope as well
%coeffyx = lregress(Y,X);
[COEF,Serr,XTXI,Rsq,Fval,TTAB,FITS,RES]= lregress(Y,X);
coeffyx = COEF;
%disp('t value is (bottom row is for slope): ')
%TTAB(:,3)
disp('coeffs., s.d. of slopes, t values (bottom row is for slope): ')
TTAB
disp('Serr is '), Serr

coeffxy = lregress(X,Y);
coeffxy = [-coeffxy(1)/coeffxy(2);1/coeffxy(2)]; 
coeff = [coeffyx coeffxy];

minx = min(X);
maxx = max(X);
t1 = linspace(minx,maxx,16);
t1 = linspace(0,maxx,16);
XI = [ones(1,16);t1];
YI = coeff'*XI;

plot(X,Y,'+',XI(2,:),YI(1,:),'linewidth',2.5)
% evd changed to just plot out regression of Y on X
%plot(X,Y,'+',XI(2,:),YI(1,:),XI(2,:),YI(2,:),':')

return
