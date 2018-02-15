clc; 
close all;

%% DFT of one slice of the image
L = 6*pi;
n = 64; 
x = linspace(0,L,n);
y = linspace(0,L,n);

[X,Y] = meshgrid(x,y);

% Z = sin(X);
Z = Y.*sin(X) - X.*cos(Y);
Zfft = fft2(Z);
linex = Z(n/2,:);
liney = Z(:,n/2); 
fftx = fft(linex);
ffty = fft(liney); 

k = 0:length(fftx)-1;
k = 1./k; % normalized wavelength

%% plotting
figure('color','white')
subplot(3,3,1)
surf(X,Y,Z)
xlabel('x')
ylabel('y')
zlabel('z')
title('Input Signal')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,2)
contourf(X,Y,abs(fftshift(Zfft)))
xlabel('x')
ylabel('y')
title('Contour of Fourier Transformed Image')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,3)
imagesc(abs(fftshift(Zfft)))
xlabel('x')
ylabel('y')
title('Contour of Fourier Transformed Image')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,4)
plot(x,linex,'linewidth',1.3,'color',rgb('MediumBlue'))
xlabel('x')
ylabel('Center line value of signal')
set(gca,'fontsize',12,'fontName','Times New Roman')


subplot(3,3,5)
scatter(k(1:n/2),abs(fftx(1:n/2)),36,'MarkerEdgeColor',rgb('MediumBlue'),...
              'MarkerFaceColor',rgb('MediumBlue'))
xlabel('Normalized Wavelength')
ylabel('Amplitude of Fourier Transform')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,6)
periodogram(linex);
set(gca,'fontsize',12,'fontName','Times New Roman')


subplot(3,3,7)
plot(y,liney,'linewidth',1.3,'color',rgb('MediumBlue'))
xlabel('y')
ylabel('Center line value of signal')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,8)
scatter(k(1:n/2),abs(ffty(1:n/2)),36,'MarkerEdgeColor',rgb('MediumBlue'),...
              'MarkerFaceColor',rgb('MediumBlue'))
xlabel('Normalized Wavelength')
ylabel('Amplitude of Fourier Transform')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,9)
periodogram(liney);
set(gca,'fontsize',12,'fontName','Times New Roman')