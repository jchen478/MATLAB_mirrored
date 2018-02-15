%%%
%%% Plot properties vs cluster size
%%%

clc;
clear;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define all relevant cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cases = {'theta0','theta1','theta3','theta6','helical'};
titleName = {'$\theta_{eq} = 0$','$\theta_{eq} = 0.1$',...
    '$\theta_{eq} = 0.3$','$\theta_{eq} = 0.6$',...
    '$(\theta_{eq},\phi_{eq})=(0.8,0.7)$'}; 
% markerArr = ['o','s','>','d','<']; 


% cases = {'theta0','theta1','theta3'}; 
% titleName = {'$\theta_{eq} = 0$','$\theta_{eq} = 0.1$','$\theta_{eq} = 0.3$'}; 
% 
% cases = {'theta6','helical'}; 
% titleName = {'$\theta_{eq} = 0.6$','$(\theta_{eq},\phi_{eq})=(0.8,0.7)$'}; 


% cases = {'theta6','helical'}; 
nCases = length(cases); 



dataFile = cell(9,1);
for j=1:9
    dataFile{j} = 'fig10_property_vs_S/'; 
    for i=1:nCases
        dataFile{j} = [dataFile{j},cases{i},'_']; 
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load all relevant cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nCases
    name = (['sub1_',cases{i},'.mat']); 
    C1(i) = load(name); 
    name = (['sub3_',cases{i},'.mat']); 
    C3(i) = load(name);
    name = (['sub7_',cases{i},'.mat']);
    C7(i) = load(name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NC vs. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{1} = [dataFile{1},'NC']; 
figure()
hold on
for i=1:nCases  
    for j=1:C1(i).nMu
        S = C3(i).S(j,:); 
        NC = C1(i).number_of_contacts(j,:); 
        scatter(S,NC,'filled');
        xlim([2 inf])        
    end
end
xlabel('$<\bar{S}>$')
ylabel('$<\bar{N}_C>$')
text(0.65, 0.3,titleName,'Units','Normalized','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I vs. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{2} = [dataFile{2},'I']; 
figure()
hold on
for i=1:nCases
    for j=1:C1(i).nMu
        S = C3(i).S(j,:); 
        intensity = C1(i).intensity(j,:); 
        scatter(S,intensity,'filled');
        xlim([2 inf])    
    end
end
xlabel('$<\bar{S}>$')
ylabel('$I$')
text(0.65, 0.3,titleName,'Units','Normalized','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sigP vs. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{3} = [dataFile{3},'sigP']; 
figure()
hold on
for i=1:nCases  
    for j=1:C1(i).nMu
        S = C3(i).S(j,:); 
        particle_stress = C1(i).particle_stress(j,:); 
        scatter(S,particle_stress,'filled');
        xlim([2 inf])        
    end
end
xlabel('$<\bar{S}>$')
ylabel('$\sigma_p L^4/ E_Y I$')
text(0.65, 0.3,titleName,'Units','Normalized','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% N1 vs. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{4} = [dataFile{4},'N1']; 
figure()
hold on
for i=1:nCases  
    for j=1:C1(i).nMu
        S = C3(i).S(j,:); 
        norm1 = C1(i).norm1(j,:); 
        scatter(S,norm1,'filled');
        xlim([2 inf])        
    end
end
xlabel('$<\bar{S}>$')
ylabel('$N_1 L^4/ E_Y I$')
text(0.65, 0.3,titleName,'Units','Normalized','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% N2 vs. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{5} = [dataFile{5},'N2']; 
figure()
hold on
for i=1:nCases  
    for j=1:C1(i).nMu
        S = C3(i).S(j,:); 
        norm2 = C1(i).norm2(j,:); 
        scatter(S,norm2,'filled');
        xlim([2 inf])        
    end
end
xlabel('$<\bar{S}>$')
ylabel('$N_2 L^4/ E_Y I$')
text(0.65, 0.3,titleName,'Units','Normalized','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eelas vs. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{6} = [dataFile{6},'Eelas']; 
figure()
hold on
for i=1:nCases  
    for j=1:C1(i).nMu
        S = C3(i).S(j,:); 
        Eelas = C1(i).Eelas(j,:); 
        scatter(S,Eelas,'filled');
        xlim([2 inf])        
    end
end
xlabel('$<\bar{S}>$')
ylabel('$E_{elas}$')
text(0.65, 0.3,titleName,'Units','Normalized','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eta_rel vs. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{7} = [dataFile{7},'etarel']; 
figure()
hold on
for i=1:nCases  
    for j=1:C1(i).nMu
        S = C3(i).S(j,:); 
        relative_viscosity = C1(i).relative_viscosity(j,:); 
        scatter(S,relative_viscosity,'filled');
        xlim([2 inf])        
    end
end
xlabel('$<\bar{S}>$')
ylabel('$\eta_{rel}$')
text(0.65, 0.3,titleName,'Units','Normalized','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dyy vs. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{8} = [dataFile{8},'Dyy']; 
figure()
hold on
for i=1:nCases  
    for j=1:C1(i).nMu
        S = C3(i).S(j,:); 
        Dyy = C7(i).Dyy(j,:); 
        scatter(S,Dyy,'filled');
        xlim([2 inf])        
    end
end
xlabel('$<\bar{S}>$')
ylabel('$D_{yy}$')
text(0.65, 0.3,titleName,'Units','Normalized','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dzz vs. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{9} = [dataFile{9},'Dzz']; 
figure()
hold on
for i=1:nCases  
    for j=1:C1(i).nMu
        S = C3(i).S(j,:); 
        Dzz = C7(i).Dzz(j,:); 
        scatter(S,Dzz,'filled');
        xlim([2 inf])        
    end
end
xlabel('$<\bar{S}>$')
ylabel('$D_{zz}$')
text(0.65, 0.3,titleName,'Units','Normalized','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:9
    figure(i)
    print(dataFile{i},'-dpng')
end