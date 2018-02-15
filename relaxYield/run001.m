%% Read files
clc; clear; close all;

ninit = 11;
ncase = 3;
startAvg = 40;

strain = linspace(1000,1500,ninit);
ind = 0:ncase-1;

filename = cell(length(strain),4,ncase);
SteadyVal = zeros(ninit,ncase,3); 
SeffArr = zeros(ncase,1);
SteadyValPreRelax = zeros(ninit,ncase,3); 

for i=1:ncase
    for j=1:ninit
        
        filename{j,1,i} = ['Parameters_case',num2str(ind(i))','.in'];
        filename{j,2,i} = ['Stress_tensor_case',num2str(ind(i))','_',num2str(strain(j)),'.txt'];
        filename{j,3,i} = ['Number_of_Contacts_case',num2str(ind(i))','_',num2str(strain(j)),'.txt'];
        filename{j,4,i} = ['Eelastic_case',num2str(ind(i))','_',num2str(strain(j)),'.txt'];
        [SteadyValPreRelax,SteadyVal,SeffArr] = relaxYield(filename{j,1,i},filename{j,2,i},filename{j,3,i},filename{j,4,i},i,j,SteadyValPreRelax,SteadyVal,startAvg,SeffArr);
    end
end

%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 20;
linewidth = 2;

%% Figure formatting
for figNum=1:ncase
    figure(figNum)
    subplot(3,1,1)
    hold on;
    ylim([0 50])
    for i=1:3
        subplot(3,1,i)
        hold on
        box on
        xlabel('\gamma')
        set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
        set(gca,'YMinorTick','on','XMinorTick','on')
        set(gca,'ticklength',tickx*get(gca,'ticklength'))
        set(gcf, 'color','white')
        xlim([-inf 40])
    end
end

SteadyValAvg = mean(SteadyVal,1)
SteadyValSig = std(SteadyVal,1)

SteadyValPreRelaxAvg = mean(SteadyValPreRelax,1);
SteadyValPreRelaxSig = std(SteadyValPreRelax,1);

%% plot steady value with 1/Seff
figure(ncase+1)
hold on
box on
errorbar(1./SeffArr,SteadyValPreRelaxAvg(:,:,1),SteadyValPreRelaxSig(:,:,1)/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',rgb('Crimson'),'Linewidth',linewidth,'color',rgb('Crimson'));
errorbar(1./SeffArr,SteadyValAvg(:,:,1),SteadyValSig(:,:,1)/2,'-.o','MarkerSize',markersize,...
    'MarkerEdgeColor',rgb('MediumBlue'),'Linewidth',linewidth,'color',rgb('MediumBlue'));
legend('Steady shear','Relaxed to zero shear','location','northwest')
xlabel('1/\it{}S_{eff}')
ylabel('\it{\sigma_{xz} L^4/ E_Y I}')
set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',tickx*get(gca,'ticklength'))
set(gcf, 'color','white')