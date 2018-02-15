close all;
clc; 

File = fopen('Sandy19_n=60_Uscale=300_mustat=0.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
g0 = fscanf(File,formatSpec,sizeA)';
fclose(File);

File = fopen('Sandy19_n=60_Uscale=300_mustat=10.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
g10 = fscanf(File,formatSpec,sizeA)';
fclose(File);

File = fopen('Sandy19_n=60_Uscale=300_mustat=100.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
g100 = fscanf(File,formatSpec,sizeA)';
fclose(File);

File = fopen('Sandy19_n=60_Uscale=300_mustat=1000.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
g1000 = fscanf(File,formatSpec,sizeA)';
fclose(File);

File = fopen('Sandy19_n=60_Uscale=300_mustat=20.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
g20 = fscanf(File,formatSpec,sizeA)';
fclose(File);

File = fopen('Sandy19_n=60_Uscale=300_mustat=200.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
g200 = fscanf(File,formatSpec,sizeA)';
fclose(File);

File = fopen('Sandy19_n=60_Uscale=300_mustat=50.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
g50 = fscanf(File,formatSpec,sizeA)';
fclose(File);

File = fopen('Sandy19_n=60_Uscale=300_mustat=500.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
g500 = fscanf(File,formatSpec,sizeA)';
fclose(File);

figure('color','white')
hold on;
plot(g0(2:end,1),g0(2:end,2),'-o')
plot(g10(2:end,1),g10(2:end,2),'-o')
plot(g100(2:end,1),g100(2:end,2),'-o')
plot(g1000(2:end,1),g1000(2:end,2),'-o')
plot(g20(2:end,1),g20(2:end,2),'-o')
plot(g200(2:end,1),g200(2:end,2),'-o')
plot(g50(2:end,1),g50(2:end,2),'-o')
plot(g500(2:end,1),g500(2:end,2),'-o')
% scatter(g1(2:end,1),g1(2:end,2))
% scatter(g5(2:end,1),g5(2:end,2))
% scatter(g10(2:end,1),g10(2:end,2))
% scatter(g20(2:end,1),g20(2:end,2))
% scatter(g50(2:end,1),g50(2:end,2))

legend('\mu = 0','\mu = 10','\mu = 100','\mu = 1000','\mu = 20'...
    ,'\mu = 200','\mu = 50','\mu = 500')
set(gca,'fontsize',18,'fontname','Times New Roman')
ylabel('\bf{g(r)}','fontsize',18)
xlabel('\bf{r/L}','fontsize',18)
title('Pair distribution (n,U_{scale}) = (60,300)','fontsize',18)
%title('N_{fib} = 160, N_{seg} = 5, \theta = 0.8, \phi = 0.7, Seff = 0.05, r_p = 75, nL^3 = 20','fontsize',12)

