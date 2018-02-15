clc; 
close all; 

%% frictionless

File = fopen('dx_frictionless.txt','r');
formatSpec = '%f';
sizeA = [Inf];
dx = fscanf(File,formatSpec,sizeA); 
fclose(File);

dx = sort(dx)
ncases = length(dx); 
S = cell(ncases,1); 
F = cell(ncases,1);

for i=1:ncases   
    S{i} = ['dx = ',num2str(dx(i),'%.2f\n')];  
    F{i} = ['Shannon_entropy_frictionless_',num2str(dx(i),'%.2f\n'),'.txt'];
end

figure('color','white')
hold on;
xlabel('strain \gamma')
ylabel('\Delta H')
set(gca,'fontsize',16)
title('Frictionless Fiber Suspension')

for i=1:ncases   
    File = fopen(F{i},'r');
    formatSpec = '%f %f';
    sizeA = [2 Inf];
    E = fscanf(File,formatSpec,sizeA)';
    fclose(File);

    time = E(2:end,1);
    entropy = E(2:end,2);
%     time = linspace(0,640,length(entropy)); 
    
    dxScale = E(1,1); 
    H_homo = E(1,2); 
    
    scatterSize = linspace(10,10,length(time));
    
    pchange = -(entropy-H_homo)/H_homo;
    scatter(time,pchange,scatterSize,'filled'); 
    
end


h = legend(S{:});
set(h,'Location','BestOutside') 


%% with friction

File = fopen('dx_friction.txt','r');
formatSpec = '%f';
sizeA = [Inf];
dx = fscanf(File,formatSpec,sizeA); 
fclose(File);

dx = sort(dx)
ncases = length(dx); 
S = cell(ncases,1); 
F = cell(ncases,1);

for i=1:ncases   
    S{i} = ['dx = ',num2str(dx(i),'%.2f\n')];  
    F{i} = ['Shannon_entropy_friction_',num2str(dx(i),'%.2f\n'),'.txt'];
end

figure('color','white')
hold on;
xlabel('strain \gamma')
ylabel('\Delta H')
set(gca,'fontsize',16)
title('Fiber Suspension with Friction')

for i=1:ncases   
    File = fopen(F{i},'r');
    formatSpec = '%f %f';
    sizeA = [2 Inf];
    E = fscanf(File,formatSpec,sizeA)';
    fclose(File);

%     time = E(2:end,1);
    entropy = E(2:end,2);
    time = linspace(0,640,length(entropy)); 
    
    dxScale = E(1,1); 
    H_homo = E(1,2); 
    
    scatterSize = linspace(10,10,length(time));
    
    pchange = -(entropy-H_homo)/H_homo;
    scatter(time,pchange,scatterSize,'filled'); 
    xlim([0 640])
    
end


h = legend(S{:});
set(h,'Location','BestOutside') 





