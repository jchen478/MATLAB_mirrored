clc; 
close all;

File = fopen('Binning_parameters.txt','r');
nbin = fscanf(File,'%d',[3 1])';
fclose(File);


File = fopen('Binning_result.txt','r');
signal = fscanf(File,formatSpec,[nbin(1),nbin(2)*nbin(3)])';
fclose(File);

signal = signal./12800;

signal = reshape(signal',[nbin(1),nbin(2), nbin(3)]); 

for i=1:nbin(3)
    
   signal(:,:,i) = signal(:,:,i)'; 
    
end

%% DFT of one slice of the image

test = signal(:,:,50); 

figure('color','white')
subplot(2,2,1)
imagesc(test);

Y = fft2(test);
subplot(2,2,2)
imagesc(abs(fftshift(Y)))


subplot(2,2,3)
periodogram(Y)

Yavg = mean(Y);
subplot(2,2,4)
periodogram(Yavg)

%% homogeneous case

File = fopen('Binning_parameters_homo.txt','r');
% formatSpec = '%d %d %d';
nbin = fscanf(File,'%d',[3 1])';
fclose(File);


File = fopen('Binning_result_homo.txt','r');
signal = fscanf(File,formatSpec,[nbin(1),nbin(2)*nbin(3)])';
fclose(File);

signal;

signal = reshape(signal',[nbin(1),nbin(2), nbin(3)]); 

for i=1:nbin(3)
    
   signal(:,:,i) = signal(:,:,i)'; 
    
end

%% DFT of one slice of the image

test = signal(:,:,50); 

figure('color','white')
subplot(2,2,1)
imagesc(test);

Y = fft2(test);
subplot(2,2,2)
imagesc(abs(fftshift(Y)))


subplot(2,2,3)
periodogram(Y)

Yavg = mean(Y);
subplot(2,2,4)
periodogram(Yavg)