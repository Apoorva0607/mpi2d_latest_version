function m2=intshft(m,sh)
[nx,ny]=size(m);
inx=(1+abs(sh(1))) : (nx-abs(sh(1)));
iny=(1+abs(sh(2))) : (ny-abs(sh(2)));
m2=m;
m2(inx-sh(1),:)= m(inx,:);      % this is original
m2(:,iny-sh(2))= m2(:,iny);
return;