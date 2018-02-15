function clusterTrajectoryFun(filename)

%% Read files
File = fopen(filename,'r');
data = fscanf(File,'%f',[2 Inf])';
fclose(File);

%% Assign data to variables
t = data(:,1);
elas = data(:,2);

%% plot
% figure()
hold on
% scatter(t,elas,20,'filled');
plot(t,elas,'-.o');

end