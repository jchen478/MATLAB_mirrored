close all;

%% declare number of files
casename = {'15', '15_run2','17','17_run2'};
legendname = {'config 1','config 2'}; 
ncase = length(casename);
filename = cell(4,1);
for i=1:ncase    
   filename{i} = ['Eelastic',casename{i},'.txt'];    
end

%% obtain dimensions and key variables
sizeA = [4 Inf];
File = fopen(filename{1},'r');
data = fscanf(File,'%f',sizeA)';
fclose(File);

strain_stop = data(1,1); 
dt = data(1,2); 
elastic_write = data(1,3);
stress_write = data(1,4);
data(1,:) = []; 

t = data(:,1);
tAfterStop = t-strain_stop; 
ndata = length(t); 
nNorm = strain_stop/elastic_write/dt;

%% allocate space for arrays
E = zeros(ndata,ncase); 
Ebend = zeros(ndata,ncase); 
Etwist = zeros(ndata,ncase); 
Enorm = zeros(ndata,ncase); 

%% read all files
for i=1:ncase
    filename{i}
    File = fopen(filename{i},'r');
    data = fscanf(File,'%f',sizeA)';
    fclose(File);
    data(1,:) = [];
    E(:,i) = data(:,4);
    Ebend(:,i) = data(:,2);
    Etwist(:,i) = data(:,3);
    Enorm(:,i) = E(:,i)./data(nNorm,4);
end

%% plot

for j=1:2:ncase
 
    % E_elastic
    figure
    subplot(2,1,1)
    hold on
    scatter(tAfterStop, E(:,j),10,'MarkerEdgeColor',rgb('Crimson'),...
        'MarkerFaceColor',rgb('Crimson'));
    scatter(tAfterStop, E(:,j+1),10,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'));
    ylabel('\it{}E_{elastic}*')
    
    % normalized elastic
%     figure
    subplot(2,1,2)
    hold on
    scatter(tAfterStop, Enorm(:,j),10,'MarkerEdgeColor',rgb('Crimson'),...
        'MarkerFaceColor',rgb('Crimson'));
    scatter(tAfterStop, Enorm(:,j+1),10,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'));
    ylabel('\it{}E_{norm}')
    
end

%% figure formatting
tickx = 1.5;

for i=1:ncase/2
    
    figure(i)
    
    for j=1:2
        subplot(2,1,j)
        hold on;
        box on;
        set(gca,'fontsize',24,'linewidth',2,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
        xlabel('\gamma')
        xlim([-20 inf])
        ylim([0 inf])
        legend(legendname{:})
    end
end