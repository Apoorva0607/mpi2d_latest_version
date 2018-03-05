function [m2,sh]=mpi_leastsquares(m,referenceFrame,testFrame,fr_num,rangex,rangey,slice)

grids=[4 3 2 1];
mframe=referenceFrame;
ww=0.5+cosw2d(size(m)); % weighting function

w3=zeros(5,5);

temp=strcat('Output/',int2str(slice),'/buffer_cinemri.mat');
load(temp);
m=[2 1 1 1];

sh=[0 0];
for n=grids
    [x,y]=ndgrid(sh(1)+n*(-2:2),sh(2)+m(n)*(-2:2));    
    for k=1:length(x(:))
        
        if(fr_num<=10)           
            try                
             %%   r=buffer_cinemri(rangey(1)+x(k):rangey(length(rangex))+x(k),rangex(1)+y(k):rangex(length(rangey))+y(k),fr_num); %%DEV TEST ON MARCH 18
                r=buffer_cinemri(:,:,fr_num);
            catch
                r=imgshft(testFrame,[x(k) y(k)]);                
            end
        else
            r=imgshft(testFrame,[x(k) y(k)]);
        end        
        r=r-mframe;        
        r=r.*ww;
        w3(k)=r(:)'*r(:);
    end    
    [~,i1]=min(w3(:));
    sh=[x(i1) y(i1)];    
end

%     disp([fr_num sh]);

m2=imgshft(testFrame,sh);


function m2=imgshft(m,sh)
m2=intshft(m,round(sh));
if(sh==0) 
    return; end;
return;


function m2=intshft(m,sh)
[nx,ny]=size(m);
inx=(1+abs(sh(1))) : (nx-abs(sh(1)));
iny=(1+abs(sh(2))) : (ny-abs(sh(2)));
m2=m;
m2(inx-sh(1),:)= m(inx,:);      % this is original
m2(:,iny-sh(2))= m2(:,iny);
return;

