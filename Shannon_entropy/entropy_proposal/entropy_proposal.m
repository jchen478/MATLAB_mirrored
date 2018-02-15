clc; 
close all;

%% with friction

ncases = 1; 
dx = [107.75];

S = cell(ncases,1); 
F = cell(ncases,1);

for i=1:ncases   
    S{i} = ['dx = ',num2str(dx(i),'%.2f\n')];  
    F{i} = ['Shannon_entropy_friction_',num2str(dx(i),'%.2f\n'),'.txt'];
end

 File = fopen(F{1},'r');
 formatSpec = '%f %f';
 sizeA = [2 Inf];
 E = fscanf(File,formatSpec,sizeA)';
 fclose(File);
 
 N = length(E(2:end,2)); 
 time640 = linspace(0,640,N); 
 
 entropy = zeros(N,4); 

for i=1:ncases   
    
    File = fopen(F{i},'r');
    formatSpec = '%f %f';
    sizeA = [2 Inf];
    E = fscanf(File,formatSpec,sizeA)';
    fclose(File);

    entropy(:,i) = E(2:end,2);

end

%% frictionless

ncases = 1; 
dx = [107.75];

S = cell(ncases,1); 
F = cell(ncases,1);

for i=1:ncases   
    S{i} = ['dx = ',num2str(dx(i),'%.2f\n')];  
    F{i} = ['Shannon_entropy_frictionless_',num2str(dx(i),'%.2f\n'),'.txt'];
end

 File = fopen(F{1},'r');
 formatSpec = '%f %f';
 sizeA = [2 Inf];
 E = fscanf(File,formatSpec,sizeA)';
 fclose(File);
 
 N = length(E(2:end,2)); 
 time640 = linspace(0,640,N); 

for i=1:ncases 
    
    File = fopen(F{i},'r');
    formatSpec = '%f %f';
    sizeA = [2 Inf];
    E = fscanf(File,formatSpec,sizeA)';
    fclose(File);

    N = length(E(2:end,2));
    
    entropy(1:N,2+i) = E(2:end,2);

end

figure('color','white')
hold on;
xlabel('\gamma')
ylabel('\it{H}')
box on
set(gca,'YMinorTick','on','XMinorTick','on')
text(340,6.25,'with friction','FontSize',22,'fontName','Times New Roman')
text(340,7.55,'without friction','FontSize',22,'fontName','Times New Roman')
set(gca,'fontsize',22,'fontName','Times New Roman','linewidth',2)
ylim([6 7.8])
xlim([0 800])
tiy=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tiy,'%.1f'))
 
scatter(time640,entropy(:,1),10,'MarkerEdgeColor',rgb('Crimson'),...
              'MarkerFaceColor',rgb('Crimson'))
scatter(time640,entropy(1:length(time640),3),10,'MarkerEdgeColor',rgb('MediumBlue'),...
              'MarkerFaceColor',rgb('MediumBlue'))
