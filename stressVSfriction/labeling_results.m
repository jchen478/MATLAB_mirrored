clc;
close all;

dataPath = '../data_stressVSfriction/clusterStrain_min2/';

markersize = 40;
linewidth = 2.5;

system = [12800];
mu = [0 1 2 3 4 5 7 10 15 17 20 23];
scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];
nmu = length(mu);
nsystem = length(system);
ncase = length(system)*length(mu);
legendArr = cell(nmu,1);

%% clusterStrain
file = cell(ncase,1);
figure()
hold on;
for j=1:nsystem
    for i=1:nmu        
        file{i+(j-1)*nmu} = [dataPath,'clusterStrain_nfib',num2str(system(j)),'_',num2str(mu(i)),'.txt'];
        File = fopen(file{i+(j-1)*nmu},'r');
        data = fscanf(File,'%f',[2 Inf])';
        fclose(File);
        legendArr{i} = ['$\mu =$ ',num2str(mu(i))];
        sort(data);        
        scatter(data(:,1),data(:,2),markersize,scatterMarker(i+(j-1)*nmu),'filled','Linewidth',linewidth)
    end
end

legend(legendArr,'location','bestoutside')
box on
xlabel('$\gamma$')
ylabel('$N_{cluster}$')
title(['$N_{fib} =$ ',num2str(system(1))])
ylim([0 inf])

