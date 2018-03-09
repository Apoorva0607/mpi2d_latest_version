for i=1:size(RegisteredData,4)
     z=fft2(squeeze(RegisteredData(80,:,:,i)));
   z1=fft2(squeeze(RegisteredData(65,:,:,i)));
   z2=fft2(squeeze(RegisteredData(75,:,:,i)));
   
    zp=[z(:,1:4) zeros(144,8) z(:,5:8)];
  zp1=[z1(:,1:4) zeros(144,8) z1(:,5:8)];
    zp2=[z1(:,1:4) zeros(144,8) z2(:,5:8)];
     G(:,:,i)=abs(ifft2(zp));
    G1(:,:,i)=abs(ifft2(zp1));
    G2(:,:,i)=abs(ifft2(zp2));
    
   
end
%%
for i=1:size(Data.Systole,4)
    figure(1);imagesc(rot90((G(:,:,i))));colormap gray;title(['Unregistered 80, time frame ',num2str(i)]);colorbar
    figure(9);imagesc(rot90((G1(:,:,i))));colormap gray;title(['Unregistered 65, time frame ',num2str(i)]);colorbar
    figure(10);imagesc(rot90((G2(:,:,i))));colormap gray;title(['Unregistered 75, time frame ',num2str(i)]);colorbar
     figure(8);imagesc(rot90(Greg(:,:,i)));colormap gray;title(['Registered ,time frame ',num2str(i)]);colorbar
    figure(3);imagesc(squeeze(RegisteredData(:,:,5,i)));colormap gray;
%     if i<size(Data.Systole,4)
% %     figure(2);imagesc(rot90(G(:,:,i+1)));colormap gray;title(['Unregistered,time frame ',num2str(i+1)]);
%     else
%        break;
%     end
    pause;
end
%%
for i=1:size(RegisteredData_new,4)
     zreg=fft2(squeeze(RegisteredData_new(75,:,:,i)));
    
    zpreg=[zreg(:,1:4) zeros(144,8) zreg(:,5:8)];
%       zpreg=[zreg(:,1:4)  zreg(:,5:8)];
    Greg(:,:,i)=abs(ifft2(zpreg));
   
end
%%
for i=1:size(S.Registered,4)
     zreg=fft2(squeeze(S.Registered(75,:,:,i)));
    
    zpreg=[zreg(:,1:4) zeros(144,8) zreg(:,5:8)];
%       zpreg=[zreg(:,1:4)  zreg(:,5:8)];
    Greg1(:,:,i)=abs(ifft2(zpreg));
   
end
%%
% Line profile for registered and unregistered
line1=squeeze(G(80,:,:));
figure(5);imagesc(line1);colormap gray;
line2=squeeze(Greg(80,:,:));
figure(6);imagesc(line2);colormap gray;