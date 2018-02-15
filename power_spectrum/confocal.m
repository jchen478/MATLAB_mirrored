function confocal(rcm,spec)

%% normalize signal
rcm = rcm./12800;

%% DFT of one slice of the image
n = length(rcm);
x = 1:n;
y = 1:n;
x = x./n*2*pi;
y = y./n*2*pi;

[X,Y] = meshgrid(x,y);

rcmfft2 = fft2(rcm);

c = 255 / log(1+max(max(rcmfft2)));
rcmfft2filtered = c*log(0.0001*rcmfft2);

%% plot
figure()
subplot(3,2,1)
contourf(X,Y,rcm)
hold on
title('Fiber centers of mass distribution')

subplot(3,2,2)
imagesc(abs(fftshift(rcmfft2filtered)))
hold on
title('Filtered Fourier transform')



subplot(3,2,5)
hold on
title('Testing line')
linex = mean(rcm,2); 
plot(1:n,linex)

subplot(3,2,6)
hold on
title('Fourier transformed test line')
Ylinex = fft(linex);
xx = 1:n;
xx = xx./n;
plot(xx(2:end/2),abs((Ylinex(2:end/2)))/n)

end