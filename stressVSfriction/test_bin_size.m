clc;
clear;
close all;

binArr = [30 40 50 64];
baseArr = [0.001465 0.003027 0.005214 0.010933];
nMu = 12;
side = 1293;
dbin = side ./binArr;

nCase = length(binArr);

binArrLegend = cell(nCase,1);

IntensityBin = zeros(nMu,4);

% read files
for i=1:nCase
    binArrLegend{i} = ['$d_{bin} = $', num2str(dbin(i))];
    name = [num2str(binArr(i)),'.txt'];
    File = fopen(name,'r');
    data = fscanf(File,'%f',[2 Inf])';
    fclose(File);
    IntensityBin(:,i) = data(:,2);
    figure(1)
    hold on
    plot(data(:,1), (data(:,2)-baseArr(i)), '-o');
    ylabel('$I-I_{base}$')
end

for j=1:1
    figure(j)
    hold on
    legend(binArrLegend,'location','best')
    xlim([0 25])
end

%% find regression line
fit = zeros(nMu,2);
% 1. Subtract base
for i=1:nCase
    IntensityBin(:,i) = IntensityBin(:,i) - baseArr(i);
end
% 2. Find dbin
dbin = side ./binArr;
% 3. Regression for every mu
for j=1:nMu
    p = polyfit(dbin',IntensityBin(j,:)',1);
    fit(j,:) = p;
end
% 4. Obtain value at rps = some value
dbinDesired = 30;
intensityDesired = zeros(nMu,1);

for j=1:nMu
    intensityDesired(j) = fit(j,2)+dbinDesired*fit(j,1);
end

% 5. Test by plotting
figure(1)
hold on
plot(data(:,1), intensityDesired, '-.o');
title(['$d_{bin} =$ ',num2str(dbinDesired)]);