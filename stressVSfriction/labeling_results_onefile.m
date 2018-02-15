clc;
clear;
close all;

system = [6400];

dataPath = '../data_stressVSfriction/meanCount_min2/';

file = [dataPath,'meanCount_nfib12800_20.txt']; 

File = fopen(file,'r');
data = fscanf(File,'%f',[4 Inf])';
fclose(File);
nstrain = length(data);
data(:,1) = [];
[sortedData, sortI] = sort(data(:,1));
sortedData = [sortedData zeros(nstrain,1) zeros(nstrain,1)];

for ii=1:nstrain
    
    sortedData(ii,2) = data(sortI(ii),2);
    sortedData(ii,3) = data(sortI(ii),3);
    
end

S = zeros(nstrain,2);
ind = 1;
for ii=1:nstrain
    S(ii,1) = sortedData(ii,1);
    if(sortedData(ii,3) ~= 0)
        
        S(ii,2) = sortedData(ii,3)/sortedData(ii,2);
        if(sortedData(ii,1) == 800)
            sortedData(ii,1)
            sortedData(ii,2)
            sortedData(ii,3)
        end
    end
end

scatter(S(:,1),S(:,2),'markerFaceColor',rgb('MediumBlue'))

box on
xlabel('$\gamma$')
ylabel('$S$')
title(['$N_{fib} =$ ',num2str(system(1))])
ylim([0 inf])