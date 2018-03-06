function [m2,shfts]=staticshft(m,grids,referenceFrame, numSkipRegistration,rangex,rangey,slice)

ww=0.5+cosw2d(size(m)); % weighting function

m2=m;
w3=zeros(5,5);

mframe=m2(:,:,referenceFrame);  % do not update reference frame after

outpath=strcat('Output/',int2str(slice));;
flname=strcat(outpath,'/buffer_cinemri.mat');
load(flname);

shfts=zeros(size(m,3),2);

m=[2 1 1 1];

for it=1+numSkipRegistration:size(m2,3)
    sh=[0 0];
    for n=grids
        
        [x,y]=ndgrid(sh(1)+n*(-2:2),sh(2)+m(n)*(-2:2));
        
        for k=1:length(x(:))
            if(it<=10)
                try
                    %r=buffer_cinemri(rangey(1)+x(k):rangey(length(rangex))+x(k),rangex(1)+y(k):rangex(length(rangey))+y(k),it);
                    r=buffer_cinemri(:,:,it); %DEV TEST ON MARCH 18
                catch
                    r=imgshft(m2(:,:,it),[x(k) y(k)]);                    
                end
            else
                r=imgshft(m2(:,:,it),[x(k) y(k)]);
            end
            
            r=r- mframe;            
            r=r.*ww;
            w3(k)=r(:)'*r(:);
        end
        
        [~,i1]=min(w3(:));
        sh=[x(i1) y(i1)];
    end    
    
%     disp([it sh]);   
    
    m2(:,:,it)=imgshft(m2(:,:,it),sh);    
    shfts(it,:)=sh;
end
temp=strcat('Output/',int2str(slice),'/shfts.mat');
save(temp,'shfts')


% shift the image on sh pixels

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
