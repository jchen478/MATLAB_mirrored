function confocal_new(Z,spec,n)

%% normalize signal
% Z = Z./12800;
% 
% for i=1:n
%     for j=1:n
%         
%         if Z(i,j) < 3
%             Z(i,j) = 0;
%         elseif Z(i,j) >= 3
%             Z(i,j) = 5;
%         end
%         
%     end
% end

%% parameters
L = length(Z); 
x = linspace(0,L,n);
y = linspace(0,L,n);
[X,Y] = meshgrid(x,y);

Zfft = fft2(Z);
% linex = Z(n/2,:);
% liney = Z(:,n/2); 
linex = mean(Z,1);
liney = mean(Z,2);
fftx = fft(linex);
ffty = fft(liney); 

k = 0:length(fftx)-1;
k = 1./k; % normalized wavelength

%% plotting
figure('color','white')
subplot(3,3,1)
[C,h] = contourf(X,Y,Z);
xlabel('x')
ylabel('y')
title(spec)
set(gca,'fontsize',12,'fontName','Times New Roman')
colormap parula 
colorbar
set(h,'LineColor','none')

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
ylabel('Averaged value of signal')
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
ylabel('Averaged value of signal')
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

end