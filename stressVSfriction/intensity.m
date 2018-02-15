clc;
clear;
close all;

dataPath = '../data_stressVSfriction/intensity/';


muArr = [0 1 2 3 4 5 7 10 15 17 20 23]; 
nfibArr = [160 1280 12800];

nMu = length(muArr); 
nNfib = length(nfibArr); 

for i=1:nNfib
    figure(i)
    hold on
    for j=1:nMu
        str = [dataPath,'nfib',num2str(nfibArr(i)),'_',num2str(muArr(j)),'_ISV.txt'];        
        File = fopen(str,'r');
        data = fscanf(File,'%f',[8 Inf])';
        fclose(File);
        
        strain = data(:,1);
        I = data(:,2); 
        
        scatter(strain,I)
        
    end
    xlabel('$\gamma$')
    ylabel('$I$')
end