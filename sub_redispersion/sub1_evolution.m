%%%% Plot the evolution of stress and number of contacts
%%%% with change in box size and volume fraction

clc;
clear;
close all;

% Define common parameters
simulation_cases;

% plotting commands
xLowLim = 0;


fig = figure('Units','Inches','Position',[1 1 3.0 3.0]);
left_color = [.5 .5 0];
right_color = rgb('MediumBlue');
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

% box
name = ['Lbox.txt'];
File = fopen(name,'r');
data = fscanf(File,'%f',[7 Inf])';
fclose(File);
boxstrain = data(:,1);
sidex = data(:,2);
sidey = data(:,3);
sidez = data(:,4);

volfrac = 1280*150*pi ./(sidex.*sidey.*sidez);
a = max(volfrac)/min(volfrac)

xUpLim = max(boxstrain);

subplot(2,1,1)
title('\boldmath{$\theta_{eq},\phi_{eq}=0.6,0,\mu = 20, a = 3$}')
hold on
% yyaxis left
plot(data(:,1),volfrac*100,'color',rgb('Black'))
xlim([xLowLim xUpLim])
ylabel('\boldmath{$\phi\ (\%)$}')
xlabel('\boldmath{$\gamma$}')
ylim([0.2 max(volfrac)*100])

subplot(2,1,2)
% yyaxis left
hold on
plot(data(:,1),volfrac*100,'color',rgb('Black'))
xlim([xLowLim xUpLim])
ylabel('\boldmath{$\phi\ (\%)$}')
xlabel('\boldmath{$\gamma$}')
ylim([0.2 max(volfrac)*100])

% number of contacts
name = ['Number_of_Contacts.txt'];
File = fopen(name,'r');
data = fscanf(File,'%f',[5 Inf])';
fclose(File);
diff = data(2,1) - data(1,1);
% make sure that strain is increasing
for ii=2:length(data)
    if (data(ii,1) < data(ii-1,1))
        data(ii,1) = data(ii-1,1) + diff;
    end
end

% block average
block = 2;
npt = floor(length(data)/block);
B = zeros(npt,1);
Bt = zeros(npt,1);
ind = 1;
for ii=1:block:length(data)-block+1
    B(ind) = mean(data(ii:ii+block-1,4));
    Bt(ind) = mean(data(ii:ii+block-1,1));
    ind = ind + 1;
end
minNC = 100;
maxNC = 0;
% find maximum and minimum
if min(B) < minNC
    minNC = min(B);
end
if max(B) > maxNC
    maxNC = max(B);
end

subplot(2,1,1)
yyaxis right
plot(Bt,B)
xlim([xLowLim xUpLim])
ylim([0 maxNC])
ylabel('$N_C$')


% stress
name = ['Stress_tensor.txt'];
File = fopen(name,'r');
data = fscanf(File,'%f',[7 Inf])';
fclose(File);
diff = data(2,1) - data(1,1);
nStep = length(data);
% make sure that strain is increasing
for ii=2:nStep
    if (data(ii,1) < data(ii-1,1))
        data(ii,1) = data(ii-1,1) + diff;
    end
end
nfib = 1280*ones(nStep,1);
nsegArr = 5*ones(nStep,1);

% remove rows from Lbox to calculate suspension conditions at config_write
ind = 1;
lboxArr = zeros(nStep,1);
data(:,1) = round(data(:,1),1);
boxstrain = round(boxstrain,1);
for ii=1:length(boxstrain)
    if boxstrain(ii) == data(ind,1)
        lboxArr(ind) = sidex(ii);
        if ind == length(data(:,1))
            break;
        end
        ind = ind+1;
    end
end

% calculate relevant parameters
[Seff,rp,nL3,L,gamma] = pCalc(nfib, nsegArr, rps, kb, lboxArr, a, EY, Imom, eta0);

% Dimensionalize stresses
[sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
    stressDim(data(:,4), zeros(nStep,1), zeros(nStep,1), ...
    zeros(nStep,1), nseg, rps, eta0, gamma, nL3);


% Non-dimensionalize
sigmap_nondim = sigmap.*L.^4/EY/Imom;

% block average
block = 5;
npt = floor(length(sigmap_nondim)/block);
B = zeros(npt,1);
Bt = zeros(npt,1);
ind = 1;
for ii=1:block:length(sigmap_nondim)-block+1
    B(ind) = mean(sigmap_nondim(ii:ii+block-1));
    Bt(ind) = mean(data(ii:ii+block-1,1));
    ind = ind + 1;
end
minStress = 100000000;
maxStress = 0;
if min(B) < minStress
    minStress = min(B);
end
if max(B) > maxStress
    maxStress = max(B);
end

subplot(2,1,2)
yyaxis right
plot(Bt,B)
xlim([xLowLim xUpLim])
ylim([minStress maxStress])
ylabel('$\sigma_{p,xz} L^4/ E_Y I$')

sub2_average_properties