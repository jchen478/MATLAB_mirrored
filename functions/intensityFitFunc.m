function [intensityDesired] = intensityFitFunc(theta,dataPath,nfib,dbinDesired,side,nMu)

% dataPath = '../data_stressVSfriction/intensityBin/';

baseFile=[dataPath,theta,'_nfib',num2str(nfib),'_baseI.txt'];
File = fopen(baseFile,'r');
data = fscanf(File,'%f',[2 Inf])';
fclose(File);
binArr = data(:,1);
% baseArr = data(:,2);
dbin = side./binArr;
nCase = length(binArr);
IntensityBin = zeros(nMu,nCase);
intensityDesired = zeros(nMu,1);
binArrLegend = cell(nCase,1);

% read files
for i=1:nCase
    binArrLegend{i} = ['$n_x = $',num2str(binArr(i)),', $d_{bin} = $', num2str(dbin(i),'%.2f')];
    name = [dataPath,theta,'_nfib',num2str(nfib),'_',num2str(binArr(i)),'I.txt'];
    File = fopen(name,'r');
    data = fscanf(File,'%f',[2 Inf])';
    fclose(File);
    IntensityBin(:,i) = data(:,2);
end

%% find regression line
fit = zeros(nMu,2);
% 1. Regression for every mu
for j=1:nMu
    p = polyfit(dbin',IntensityBin(j,:),1);
    fit(j,:) = p;
end
% 3. Obtain value at rps = some value
for j=1:nMu
    intensityDesired(j) = fit(j,2)+dbinDesired*fit(j,1);
end
% 4. Test by plotting
% figure(figInd)
% hold on
% plot(data(:,1), intensityDesired, '-.o');
% title(['$N_{fib} = $ ',num2str(nfib),' $d_{fit} =$ ',num2str(dbinDesired,'%.2f')]);
% title(['$N_{fib} = $ ',num2str(nfib)]);
end