function [ period, power ] = stress_signal_process(x,topk)

% sample 5 points over unit strain
fs = 5;

% shift to zero-mean
x = x - mean(x);

% number of samples (only even number)
n = length(x);
if (mod(n,2) == 1)
    x(end) = []; 
    n = n-1; 
end   

y = fft(x);
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% P2 = abs(y/n); 

% figure()
% box on
% plot(f,power)
% xlabel('Frequency')
% ylabel('Power')

y0 = fftshift(y);         % shift y values
f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
power0 = abs(y0/n).^2;    % 0-centered power

% extract only positive frequency range
f0(1:n/2) = []; 
power0(1:n/2) = [];

% figure('Units','Inches','Position',[1 1 4.25 4.25])
% box on
% plot(f0,power0)
% xlabel('Frequency')
% ylabel('Power')

[B,I] = sort(power0,'descend');

power = B(1:topk); 
period = zeros(topk,1);

for k=1:topk
    fk = f0(I(k));
    period(k) = 1/fk;
end

end