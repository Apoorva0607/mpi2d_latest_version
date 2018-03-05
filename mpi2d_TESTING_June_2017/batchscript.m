clc
close all
clear all

folder = '/v/raid1a/dlikhite/MRIdata/Cardiac/Verio/P012712/Processing';
%folder = '/home/mirl/bmatthew/Desktop/test/Files you should get/Processing';
%8. Run this script

cd(folder);
if(exist('mpi2d.par')) delete('mpi2d.par'); end
files = dir([folder '/*.par']);
for i=1:length(files)
    copyfile(char(files(i).name),'mpi2d.par');
    disp(['i = ' num2str(i)]);
    disp(['Processing ' files(i).name]);
    %2.6 cross-correlation registration with gaussian centering
        %if this registration screws up, run mpi2d stg=2.11 to apply zero
        %shifts then do 3.11 to do manual registration
    %3.12 automatic segmentation
    %3.2 extract curves
    %3.4 convert curves to deltaSIcurves
    %4.1 compute 2-compartment model params
    %6.1 display k-trans
      try
    mpi2d stg=[9.2]
    
    catch
        continue
    end
end

clear all
close all
clc

folder = '/v/raid1a/dlikhite/MRIdata/Cardiac/Verio/P021712/Processing';
%folder = '/home/mirl/bmatthew/Desktop/test/Files you should get/Processing';
%8. Run this script

cd(folder);
if(exist('mpi2d.par')) delete('mpi2d.par'); end
files = dir([folder '/*.par']);
for i=1:length(files)
    copyfile(char(files(i).name),'mpi2d.par');
    disp(['i = ' num2str(i)]);
    disp(['Processing ' files(i).name]);
    %2.6 cross-correlation registration with gaussian centering
        %if this registration screws up, run mpi2d stg=2.11 to apply zero
        %shifts then do 3.11 to do manual registration
    %3.12 automatic segmentation
    %3.2 extract curves
    %3.4 convert curves to deltaSIcurves
    %4.1 compute 2-compartment model params
    %6.1 display k-trans
      try
    mpi2d stg=[9.2]
    
    catch
        continue
    end
end

clear all
close all
clc

mail