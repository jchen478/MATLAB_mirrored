%%%
%%% Find the slopes of all properties w.r.t. system size
%%% -- should run sub1_averaged_values prior to this script
%%%

clc;
clear;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define all relevant cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cases = {'theta0','theta1','theta3','theta6','helical'}; 
cases = {'theta0','theta1','theta3'}; 
nCases = length(cases); 

dataFile = cell(9,1);
for j=1:9
    dataFile{j} = 'fig9_property_slope/'; 
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
    name = (['sub7_',cases{i},'.mat']);
    C7(i) = load(name); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slope for NC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{1} = [dataFile{1},'NC']; 
figure()
hold on
for i=1:nCases    
    slope = zeros(C1(i).nMu,1); 
    for j=1:C1(i).nMu
        slope(j) = regress((C1(i).number_of_contacts(j,:))',(C1(i).lboxArr)');
    end
    plot(C1(i).muArr,slope,'-.o')
end
legend(cases,'location','best');
xlabel('$\mu$')
ylabel('$dN_C/dL_{box}$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slope for I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{2} = [dataFile{2},'I']; 
figure()
hold on
for i=1:nCases    
    slope = zeros(C1(i).nMu,1); 
    for j=1:C1(i).nMu
        slope(j) = regress((C1(i).intensity(j,:))',(C1(i).lboxArr)');
    end
    plot(C1(i).muArr,slope,'-.o')
end
legend(cases,'location','best');
xlabel('$\mu$')
ylabel('$dI/dL_{box}$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slope for sigP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{3} = [dataFile{3},'sigP']; 
figure()
hold on
for i=1:nCases    
    slope = zeros(C1(i).nMu,1); 
    for j=1:C1(i).nMu
        slope(j) = regress((C1(i).particle_stress(j,:))',(C1(i).lboxArr)');
    end
    plot(C1(i).muArr,slope.*C1(i).L.^4/C1(i).EY/C1(i).Imom/2,'-.o')
end
legend(cases,'location','best');
xlabel('$\mu$')
ylabel('$d\sigma_p/dL_{box} * L^4/ E_Y I$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slope for N1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{4} = [dataFile{4},'N1'];
figure()
hold on
for i=1:nCases    
    slope = zeros(C1(i).nMu,1); 
    for j=1:C1(i).nMu
        slope(j) = regress((C1(i).norm1(j,:))',(C1(i).lboxArr)');
    end
    plot(C1(i).muArr,slope.*C1(i).L.^4/C1(i).EY/C1(i).Imom/2,'-.o')
end
legend(cases,'location','best');
xlabel('$\mu$')
ylabel('$dN_1/dL_{box} * L^4/ E_Y I$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slope for N2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{5} = [dataFile{5},'N2']; 
figure()
hold on
for i=1:nCases    
    slope = zeros(C1(i).nMu,1); 
    for j=1:C1(i).nMu
        slope(j) = regress((C1(i).norm2(j,:))',(C1(i).lboxArr)');
    end
    plot(C1(i).muArr,slope.*C1(i).L.^4/C1(i).EY/C1(i).Imom/2,'-.o')
end
legend(cases,'location','best');
xlabel('$\mu$')
ylabel('$dN_2/dL_{box} * L^4/ E_Y I$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slope for Eelas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{6} = [dataFile{6},'Eelas']; 
figure()
hold on
for i=1:nCases    
    slope = zeros(C1(i).nMu,1); 
    for j=1:C1(i).nMu
        slope(j) = regress((C1(i).Eelas(j,:))',(C1(i).lboxArr)');
    end
    plot(C1(i).muArr,slope,'-.o')
end
legend(cases,'location','best');
xlabel('$\mu$')
ylabel('$dE_{elas}/dL_{box}$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slope for Eta_rel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{7} = [dataFile{7},'etarel']; 
figure()
hold on
for i=1:nCases    
    slope = zeros(C1(i).nMu,1); 
    for j=1:C1(i).nMu
        slope(j) = regress((C1(i).relative_viscosity(j,:))',(C1(i).lboxArr)');
    end
    plot(C1(i).muArr,slope,'-.o')
end
legend(cases,'location','best');
xlabel('$\mu$')
ylabel('$d\eta_{rel}/dL_{box}$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slope for Dyy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{8} = [dataFile{8},'Dyy']; 
figure()
hold on
for i=1:nCases    
    slope = zeros(C7(i).nMu,1); 
    for j=1:C7(i).nMu
        slope(j) = regress((C7(i).Dyy(j,:))',(C7(i).lboxArr)');
    end
    plot(C7(i).muArr,slope,'-.o')
end
legend(cases,'location','best');
xlabel('$\mu$')
ylabel('$dD_{yy}/dL_{box}$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slope for Dzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile{9} = [dataFile{9},'Dzz']; 
figure()
hold on
for i=1:nCases    
    slope = zeros(C7(i).nMu,1); 
    for j=1:C7(i).nMu
        slope(j) = regress((C7(i).Dzz(j,:))',(C7(i).lboxArr)');
    end
    plot(C7(i).muArr,slope,'-.o')
end
legend(cases,'location','best');
xlabel('$\mu$')
ylabel('$dD_{zz}/dL_{box}$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:9
    figure(i)
    print(dataFile{i},'-dpng')
end