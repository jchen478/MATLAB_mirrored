clc; 
close all;

%% parameters
File = fopen('Binning_parameters.txt','r');
nbin = fscanf(File,'%d',[3 1])';
fclose(File);

%% xy view
File = fopen('Binning_result_xy.txt','r');
signal = fscanf(File,'%f',[nbin(1),nbin(2)])';
fclose(File);

confocal_new(signal,'xy view flocculated',nbin(1))

%% xz view
File = fopen('Binning_result_xz.txt','r');
signal = fscanf(File,'%f',[nbin(1),nbin(3)])';
fclose(File);

confocal_new(signal,'xz view flocculated',nbin(1))

%% yz view
File = fopen('Binning_result_yz.txt','r');
signal = fscanf(File,'%f',[nbin(2),nbin(3)])';
fclose(File);

confocal_new(signal,'yz view flocculated',nbin(1))

%% xy view
File = fopen('Binning_result_homo_xy.txt','r');
signal = fscanf(File,'%f',[nbin(1),nbin(2)])';
fclose(File);

confocal_new(signal,'xy view homogeneous',nbin(1))
% 
% %% xz view
% File = fopen('Binning_result_homo_xz.txt','r');
% signal = fscanf(File,'%f',[nbin(1),nbin(3)])';
% fclose(File);
% 
% confocal_new(signal,'xz view homogeneous')
% %% yz view
% File = fopen('Binning_result_homo_yz.txt','r');
% signal = fscanf(File,'%f',[nbin(2),nbin(3)])';
% fclose(File);
% 
% confocal_new(signal,'yz view homogeneous')
