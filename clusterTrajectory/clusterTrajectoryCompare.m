close all
clc
clear

id = [2 11 15 19 38 40 48]; 
% id = [2 11 15]; 

File = fopen('Cluster_results.txt','r');
nfibC = fscanf(File,'%f',[2 max(id)+1])';
fclose(File);
nfibC = nfibC(:,2); 

ncase = length(id); 
legendArr = cell(1);
filename = cell(1,ncase);

temp = zeros(ncase,2);

temp(:,2) = id;
temp(:,1) = nfibC(id(:)+1);
temp = sortrows(temp);
id = temp(:,2);

figure()
for k=1:ncase
    
    filename{k} = ['flocElastic',num2str(id(k)),'.txt'];
    clusterTrajectoryFun(filename{k}); 
    legendArr{k} = ['Final cluster size = ',num2str(nfibC(id(k)+1))];
    
end

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 24;
linewidth = 3;

for i=1:1
    figure(i)
    hold on
    box on
    set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
    set(gca,'YMinorTick','on','XMinorTick','on')
    set(gca,'ticklength',tickx*get(gca,'ticklength'))
    set(gcf, 'color','white')
    xlim([0 inf])
    xlabel('\gamma')
    ylabel('Averaged \it{}E_{elastic}')
    legend(legendArr{:},'location','best')
    xlim([0.5 inf])
end
