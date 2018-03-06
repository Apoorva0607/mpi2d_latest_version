clear b
figure(1); clf;

load ../sos_all_MID_581.mat
MID=581;

all_dia=zeros(size(sos_all));
all_sys=zeros(size(sos_all));

sos_all=sos_all(21:108,21:108,:,:);

sl=squeeze(sos_all(:,:,:,1));

sos_all=sos_all/2;


for islice=1:5
    filename=strcat('new_dia_',int2str(islice),'_MID_',int2str(MID),'.mat');
    load(filename)
    all_dia(:,:,:,islice)=new_dia_sl;
end

all_dia=abs(all_dia)/8;
all_dia=all_dia(21:108,21:108,:,:);


for islice=1:5
    filename=strcat('new_sys_',int2str(islice),'_MID_',int2str(MID),'.mat');
    load(filename)
    all_sys(:,:,:,islice)=new_sys_sl;
end

all_sys=abs(all_sys)/8;
all_sys=all_sys(21:108,21:108,:,:);

b(:,:,1,1)=sl(:,:,10);
b(:,:,1,2)=sl(:,:,10);
b(:,:,1,3)=sl(:,:,10);
b(:,:,1,4)=sl(:,:,10);
b(:,:,1,5)=sl(:,:,10);

b(:,:,1,6)=all_sys(:,:,10);
b(:,:,1,7)=all_sys(:,:,10);
b(:,:,1,8)=all_sys(:,:,10);
b(:,:,1,9)=all_sys(:,:,10);
b(:,:,1,10)=all_sys(:,:,10);

% b(:,:,1,11)=all_sys(:,:,10);
% b(:,:,1,12)=all_sys(:,:,10);
% b(:,:,1,13)=all_sys(:,:,10);
% b(:,:,1,14)=all_sys(:,:,10);
% b(:,:,1,15)=all_sys(:,:,10);

z=montage(b/2.8,'size',[2 5]);brighten(0.2)
% montage can also be used

set(gcf, 'Color' ,'w'); % set background to be white
f = getframe(figure(1)); % you get the whole picture including the axis
[im,map] = rgb2ind(f.cdata,1024);

% [im,map] = rgb2ind(f.cdata,1024,'nodither');
% [im,map] = gray2ind(f.cdata,512);

offset=10+10+8;

for k = 1:80
    k
        
    b(:,:,1,1)=sos_all(:,:,k+offset,3);
    b(:,:,1,2)=sos_all(:,:,k+offset,5)/1.35;
    b(:,:,1,3)=sos_all(:,:,k+offset,2);
    b(:,:,1,4)=sos_all(:,:,k+offset,4)/1.35;
    b(:,:,1,5)=sos_all(:,:,k+offset,1)*1.5;
    
    
    b(:,:,1,6)=all_sys(:,:,k+offset,3);    
    b(:,:,1,7)=all_sys(:,:,k+offset,5)/1.35;    
    b(:,:,1,8)=all_sys(:,:,k+offset,2);    
    b(:,:,1,9)=all_sys(:,:,k+offset,4)/1.35;    
    b(:,:,1,10)=all_sys(:,:,k+offset,1);
    
    
%     b(:,:,1,11)=all_sys(:,:,k+offset,3);    
%     b(:,:,1,12)=all_sys(:,:,k+offset,5);    
%     b(:,:,1,13)=all_sys(:,:,k+offset,2);    
%     b(:,:,1,14)=all_sys(:,:,k+offset,4);    
%     b(:,:,1,15)=all_sys(:,:,k+offset,1);       
    
    
    montage(b/9.8,'size',[2 5]);brighten(0.2)
    
    set(gcf, 'Color' ,'w');
    axis tight
    set(gca,'nextplot','replacechildren','visible','off')
    f = getframe(figure(1)); % you get the whole picture including the axis
    im(:,:,1,k) = rgb2ind(f.cdata,map);
    
    %   im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
    %   im(:,:,1,k) = gray2ind(f.cdata,map);
    pause(0.1)
end

imwrite(im,map,'Test.gif','DelayTime',0.1)

