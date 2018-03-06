clear

try
    matlabpool local 12
end

% load ../sos_all_MID_581.mat
% MID=581;
% series_start=19;

load ../sos_all_MID_603.mat
MID=603;
series_start=39;

nof=size(sos_all,3);
orig_times=(1:size(sos_all,3));

for islice=1:5
    
    sl=squeeze(sos_all(:,:,:,islice));    
    
    % diastole
        
    filename=strcat('/v/raid1/ed/MRIdata/Cardiac/Verio/P082211/ReconData/mat_files/listofselectedFrames_roi_sum_dd_Neardiastolic_series',int2str(series_start+islice-1),'000slice',int2str(islice),'.mat');
    load(filename)
%     peakloc=peakloc+1;
            
    filename=strcat('registered_dia_',int2str(islice),'_MID_',int2str(MID),'.mat');
    load(filename)
    dia_sl=registered;
    
    [sx sy sz]=size(dia_sl);
    
    dia_curves=reshape(dia_sl,[sx*sy sz]);
    dia_times=peakloc;
    dia_times=dia_times(6:end);
    
    new_dia_curves=zeros(sx*sy,nof);
    
    parfor i=1:size(dia_curves,1)
        new_dia_curves(i,:)=interp1(dia_times,dia_curves(i,:),orig_times,'cubic');
    end
    
    new_dia_sl=reshape(new_dia_curves,[sx sy nof]);    
    
    for i=1:nof
        subplot(1,2,1),imagesc(fliplr(rot90(sl(:,:,i),-1))),colormap gray,brighten(0.3)
        subplot(1,2,2),imagesc(new_dia_sl(:,:,i)),colormap gray,brighten(0.3)
        pause(0.1)
    end
    
    filename=strcat('new_dia_',int2str(islice),'_MID_',int2str(MID),'.mat');
    save(filename,'new_dia_sl')    
        
    % systole
        
    filename=strcat('/v/raid1/ed/MRIdata/Cardiac/Verio/P082211/ReconData/mat_files/listofselectedFrames_roi_sum_dd_Nearsystolic_series',int2str(series_start+islice-1),'000slice',int2str(islice),'.mat');
    load(filename)
%     peakloc=peakloc+1;
       
    filename=strcat('registered_sys_',int2str(islice),'_MID_',int2str(MID),'.mat');
    load(filename)
    sys_sl=registered;
    
    [sx sy sz]=size(sys_sl);
    
    sys_curves=reshape(sys_sl,[sx*sy sz]);
    sys_times=peakloc;
    sys_times=sys_times(6:end);

    new_sys_curves=zeros(sx*sy,nof);
    
    parfor i=1:size(sys_curves,1)
        new_sys_curves(i,:)=interp1(sys_times,sys_curves(i,:),orig_times,'cubic');
    end
    
    new_sys_sl=reshape(new_sys_curves,[sx sy nof]);    
    
    for i=1:nof
        subplot(1,2,1),imagesc(flipud(rot90(sl(:,:,i)))),colormap gray,brighten(0.3)
        subplot(1,2,2),imagesc(new_sys_sl(:,:,i)),colormap gray,brighten(0.3)
        pause(0.1)
    end
        
    filename=strcat('new_sys_',int2str(islice),'_MID_',int2str(MID),'.mat');
    save(filename,'new_sys_sl')
    
end

