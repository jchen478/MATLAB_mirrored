clc;
close all;
clear;

dataPath = '../data_stressVSfriction/meanCount_min2/'; 

markersize = 50;
tickx = 1.5;
fontsize = 24;
linewidth = 2.5;

system = [12800];
mu = [0 1 2 3 4 5 7 10 15 17 20 23];


scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];

nmu = length(mu);
nsystem = length(system);
ncase = length(system)*length(mu);
legendArr = cell(nmu,1);

CountAverage = zeros(nmu,3);
nAvg = 40;

file = cell(ncase,1);
figure()
hold on;
for j=1:nsystem
    for i=1:nmu
        
        file{i+(j-1)*nmu} = [dataPath,'meanCount_nfib',num2str(system(j)),'_',num2str(mu(i)),'.txt'];
        File = fopen(file{i+(j-1)*nmu},'r');
        dataUnsorted = fscanf(File,'%f',[4 Inf])';
        fclose(File);
        
        %% remove friction coeffient column
        dataUnsorted(:,1) = [];
        
        %% number of sample points
        nstrain = length(dataUnsorted(:,1));
        
        %% sort data
        data = zeros(size(dataUnsorted));
        [data(:,1), sortI] = sort(dataUnsorted(:,1));
        for ii=1:nstrain
            data(ii,2) = dataUnsorted(sortI(ii),2);
            data(ii,3) = dataUnsorted(sortI(ii),3);
        end
        
        %% calculate mean size S
        S = zeros(nstrain,2);
        S(:,1) = data(:,1);
        for ii=1:nstrain
            if(data(ii,2) ~= 0)
                S(ii,2) = data(ii,3)/data(ii,2);
            end
        end
        
        %% calculate average S
        CountAverage(i,1) = mu(i);
        CountAverage(i,2) = mean(S(end-nAvg:end));  % average
        CountAverage(i,3) = std(S(end-nAvg:end))/sqrt(nAvg); % standard deviation of mean
        
        %% define legend for plotting
        legendArr{i} = ['$\mu =$ ',num2str(mu(i))];
        
        %% plot S evolution with strain
        figure(1)
        hold on
        scatter(S(:,1),S(:,2),markersize,scatterMarker(i),'filled')
        figure(2)
        hold on
        scatter(S(:,1),S(:,2)/S(1,2),markersize,scatterMarker(i),'filled')
        
    end
end

% str=[num2str(system(1)),'_meanClusterSize.mat'];
% save(str,'CountAverage')

figure(1)
box on
legend(legendArr,'location','bestoutside')
xlabel('$\gamma$')
ylabel('$S$')
title(['$N_{fib} =$ ',num2str(system(1))])
ylim([0 inf])
xlim([0 inf])
set(gca,'yscale','log')

figure(2)
legend(legendArr,'location','bestoutside')
xlabel('$\gamma$')
ylabel('$S / S(\gamma = 0)$')
title(['$N_{fib} =$ ',num2str(system(1))])
ylim([0 inf])
xlim([0 inf])
% set(gca,'yscale','log')

figure(3)
errorbar(CountAverage(:,1),CountAverage(:,2),...
    CountAverage(:,3)/2,'-.o','MarkerSize',10,...
    'MarkerEdgeColor',rgb('MediumBlue'),...
    'Linewidth',2.5,'color',rgb('MediumBlue'));
xlabel('$\mu$')
ylabel('Average $S$')
title(['$N_{fib} =$ ',num2str(system(1))])
ylim([0 inf])
xlim([0 25])

%% intensity
file = cell(ncase,1);
intensityAverage = zeros(nmu,3);

dataPath = '../data_stressVSfriction/intensity/'; 

figure(4)
hold on;
for j=1:nsystem
    for i=1:nmu
        file{i+(j-1)*nmu} = [dataPath,'nfib',num2str(system(j)),'_',num2str(mu(i)),'_ISV.txt'];
        File = fopen(file{i+(j-1)*nmu},'r');
        data = fscanf(File,'%f',[8 Inf])';
        fclose(File);
        time = data(:,1);
        intensity = data(:,2);
        legendArr{i} = ['$\mu =$ ',num2str(mu(i))];
        scatter(time,intensity,7,scatterMarker(i+(j-1)*nmu),'filled','Linewidth',linewidth)
        %% calculate average I
        intensityAverage(i,1) = mu(i);
        intensityAverage(i,2) = mean(intensity(end-nAvg:end));  % average
        intensityAverage(i,3) = std(intensity(end-nAvg:end))/sqrt(nAvg); % standard deviation of mean
    end
end
legend(legendArr,'location','bestoutside')
box on
xlabel('$\gamma$')
ylabel('$I$')
set(gca,'yscale','log')
title(['$N_{fib} =$ ',num2str(system(1))])

figure(5)
errorbar(intensityAverage(:,1),intensityAverage(:,2),...
    intensityAverage(:,3)/2,'-.o','MarkerSize',5,...
    'MarkerEdgeColor',rgb('MediumBlue'),...
    'Linewidth',2.5,'color',rgb('MediumBlue'));
xlabel('$\mu$')
ylabel('Average $I$')
title(['$N_{fib} =$ ',num2str(system(1))])
ylim([0 inf])
xlim([0 25])

%% Compare S with I
Snorm = CountAverage(:,2)-min(CountAverage(:,2));
se_Snorm = CountAverage(:,3) / max(Snorm);
Snorm = Snorm / max(Snorm);
Inorm = intensityAverage(:,2)-min(intensityAverage(:,2));
se_Inorm = intensityAverage(:,3) / max(Inorm);
Inorm = Inorm / max(Inorm);

figure(6)
subplot(2,1,1)
hold on
errorbar(mu,Snorm,se_Snorm/2,...
    '-.o','MarkerSize',10,...
    'MarkerEdgeColor',rgb('Crimson'),...
    'Linewidth',2.5,'color',rgb('Crimson'));
errorbar(mu,Inorm,se_Inorm/2,...
    '-.o','MarkerSize',10,...
    'MarkerEdgeColor',rgb('MediumBlue'),...
    'Linewidth',2.5,'color',rgb('MediumBlue'));
xlabel('$\mu$')
title(['$N_{fib} =$ ',num2str(system(1))])
ylim([0 inf])
xlim([0 25])
legend('scaled $\bar{S}$','scaled $\bar{I}$','location','best')


subplot(2,1,2)
hold on
plot(mu,Snorm./Inorm,...
    '-.o','MarkerSize',10,...
    'MarkerEdgeColor',rgb('MediumBlue'),...
    'Linewidth',2.5,'color',rgb('MediumBlue'));
xlabel('$\mu$')
ylabel('scaled $\bar{S}/\bar{I}$')
ylim([0 inf])
xlim([0 25])


spreadfigures();