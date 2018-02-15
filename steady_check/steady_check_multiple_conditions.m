clc;
clear;
close all;

%% common parameters
nfibArr = [12800 ];
muArr = [0 1 2 3 4 5 7 10 15 17 20 23];
thetaArr = [1];
rpFiber = 75;
scatterMarker=['o' 'v' '>' '<' '^' 's' 'd' 'v' '>' '<' '^' 's' 'd' 'o' ];
colorArr = {rgb('DarkRed') rgb('Crimson') rgb('OrangeRed') rgb('Orange') rgb('Gold')...
    rgb('Lime') rgb('Olive') ...
    rgb('DarkGreen') rgb('LightSkyBlue') ...
    rgb('MediumBlue')...
    rgb('Plum') rgb('Purple') };

nTheta = length(thetaArr);
nMu = length(muArr);
nNfib = length(nfibArr);
muLegendArr = cell(nMu,1);
thetaNfibLegendArr = cell(nTheta*nNfib,1);
markersize = 50;

for i=1:nMu
    muLegendArr{i} = ['$\mu =$ ',num2str(muArr(i))];
end
for i=1:nTheta
    for j=1:nNfib
        thetaNfibLegendArr{(i-1)*nNfib+j} = ['$(\theta_{eq},N_{fib}) =$ (0.',num2str(thetaArr(i)),', ',num2str(nfibArr(j)),')'];
    end
end

figStart = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intensity and number of contacts vs. strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = '../data_stressVSfriction/intensity/';

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.5 0.6])
        hold on;
        ylabel('$I$')
        title(['$(\theta_{eq},N_{fib}) = (0.$',num2str(thetaArr(k)),', ',num2str(nfibArr(j)),'$)$'])
%         for i=1:nMu
%             % read files
%             name = [dataPath,'theta',num2str(thetaArr(k)),'_ISV_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
%             File = fopen(name,'r');
%             data = fscanf(File,'%f',[8 Inf])';
%             fclose(File);
%             % plot intensity evolution with strain
%             scatter(data(:,1),data(:,2),markersize,scatterMarker(i),...
%                 'MarkerFaceColor',colorArr{i},...
%                 'MarkerEdgeColor',colorArr{i})
%         end
    end
end

dataPath = '../data_stressVSfriction/nc/';

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.5 0.6])
        hold on;
        ylabel('$N_c$')
        title(['$(\theta_{eq},N_{fib}) = (0.$',num2str(thetaArr(k)),', ',num2str(nfibArr(j)),'$)$'])
        for i=1:nMu
            % read files
            name = [dataPath,'theta',num2str(thetaArr(k)),'_NC_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            data = fscanf(File,'%f',[5 Inf])';
            fclose(File);
            diff = data(2,1) - data(1,1);
            % make sure that strain is increasing
            for ii=2:length(data)
                if (data(ii,1) < data(ii-1,1))
                    data(ii,1) = data(ii-1,1) + diff;
                end
            end
            % plot number of contacts evolution with strain
            scatter(data(:,1),data(:,4),markersize,scatterMarker(i),...
                'MarkerFaceColor',colorArr{i},...
                'MarkerEdgeColor',colorArr{i})
        end
    end
end

dataPath = '../data_stressVSfriction/stress/';

for k=1:nTheta
    for j=1:nNfib
        figure('units','normalized','outerposition',[0.1 0.1 0.5 0.6])
        hold on;
        ylabel('Stress output')
        title(['$(\theta_{eq},N_{fib}) = (0.$',num2str(thetaArr(k)),', ',num2str(nfibArr(j)),'$)$'])
        for i=1:nMu
            % read files
            name = [dataPath,'theta',num2str(thetaArr(k)),'_stress_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            File = fopen(name,'r');
            data = fscanf(File,'%f',[7 Inf])';
            fclose(File);
            diff = data(2,1) - data(1,1);
            % make sure that strain is increasing
            for ii=2:length(data)
                if (data(ii,1) < data(ii-1,1))
                    data(ii,1) = data(ii-1,1) + diff;
                end
            end
            % plot number of contacts evolution with strain
            scatter(data(:,1),data(:,4),markersize,scatterMarker(i),...
                'MarkerFaceColor',colorArr{i},...
                'MarkerEdgeColor',colorArr{i})
        end
    end
end

for i=figStart:3*figStart+nNfib*nTheta-1
    figure(i)
    box on
    xlabel('$\gamma$')
    xlim([0 inf])
    legend(muLegendArr,'location','bestoutside')
    name = ['figures/steady_check_figure',num2str(i)];
    savefig([name,'.fig']);
    print(name,'-dpng')
end

% spreadfigures();
figStart = figStart + 3*nNfib*nTheta;
