clc; 
close all;

%% sine functions
L = 1; 
n = 1024
x = linspace(0,L,n); 

freq1 = 1/L; 
omega1 = 2*pi*freq1;
f1 = sin(omega1*x);
f1fft = fft(f1);

freq2 = 2/L;
omega2 = 2*pi*freq2;
f2 = sin(omega2*x);
f2fft = fft(f2); 

freq3 = 4/L;
omega3 = 2*pi*freq3;
f3 = sin(omega3*x);
f3fft = fft(f3); 

k = 0:length(f2fft)-1;
k = 1./k; % normalized wavelength


%% plot
figure('color','white')
subplot(3,3,1)
plot(x,f1,'linewidth',1.3,'color',rgb('MediumBlue'))
xlabel('Normalized Distance')
ylabel('Signal')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,2)
scatter(k(1:n/2),abs(f1fft(1:n/2)),36,'MarkerEdgeColor',rgb('MediumBlue'),...
              'MarkerFaceColor',rgb('MediumBlue'))
xlabel('Normalized Wavelength')
ylabel('Amplitude of Fourier Transform')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,3)
periodogram(f1);
set(gca,'fontsize',12,'fontName','Times New Roman')


subplot(3,3,4)
plot(x,f2,'linewidth',1.3,'color',rgb('MediumBlue'))
xlabel('Normalized Distance')
ylabel('Signal')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,5)
scatter(k(1:n/2),abs(f2fft(1:n/2)),36,'MarkerEdgeColor',rgb('MediumBlue'),...
              'MarkerFaceColor',rgb('MediumBlue'))
xlabel('Normalized Wavelength')
ylabel('Amplitude of Fourier Transform')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,6)
periodogram(f2);
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,7)
plot(x,f3,'linewidth',1.3,'color',rgb('MediumBlue'))
xlabel('Normalized Distance')
ylabel('Signal')
set(gca,'fontsize',12,'fontName','Times New Roman')


subplot(3,3,8)
scatter(k(1:n/2),abs(f3fft(1:n/2)),36,'MarkerEdgeColor',rgb('MediumBlue'),...
              'MarkerFaceColor',rgb('MediumBlue'))
xlabel('Normalized Wavelength')
ylabel('Amplitude of Fourier Transform')
set(gca,'fontsize',12,'fontName','Times New Roman')

subplot(3,3,9)
periodogram(f3);
set(gca,'fontsize',12,'fontName','Times New Roman')

