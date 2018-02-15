function [intensityDesired] = intensityFitFunc(theta,nfib,dbinDesired,side,nMu,figInd)

dataPath = '../data_stressVSfriction/intensityBin/';

baseFile=[dataPath,theta,'_nfib',num2str(nfib),'_baseI.txt'];
File = fopen(baseFile,'r');
data = fscanf(File,'%f',[2 Inf])';
fclose(File);
binArr = data(:,1);
baseArr = data(:,2);

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
%     data
    IntensityBin(:,i) = data(:,2);
%     figure(figInd)
%     hold on
%     plot(data(:,1), (data(:,2)), '-o');
end

% 
% figure(figInd)
% hold on
% xlabel('$\mu$')
% ylabel('$I$')
% legend(binArrLegend,'location','bestoutside')
% xlim([0 25])


%% find regression line
fit = zeros(nMu,2);
% 1. Subtract base
% for i=1:nCase
%     IntensityBin(:,i) = IntensityBin(:,i) - baseArr(i);
% end
% 2. Regression for every mu
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