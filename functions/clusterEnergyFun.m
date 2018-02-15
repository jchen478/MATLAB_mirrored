function clusterEnergyFun(filename)

%% Read files
File = fopen(filename,'r');
data = fscanf(File,'%f',[3 Inf])';
fclose(File);

%% Assign data to variables
% id = data(:,1);
nfibC = data(:,2);
elas = data(:,3);

%% plot
% close all
% figure()
hold on
scatter(nfibC,elas,50,'filled');

end