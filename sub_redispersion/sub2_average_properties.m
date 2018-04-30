%% Determine stress averaged over specified intervals

%%%% Plot the evolution of stress and number of contacts
%%%% with change in box size and volume fraction

clc;
clear;
close all;

% Define common parameters
% simulation_cases;
a = 1.60E-05;       % radius (m)
Imom = pi*a^4/4;    % area moment (m^4)
EY = 1e9;           % Young's modulus (Pa m^2)
eta0 = 1;           % fluid viscosity (Pa s)
nseg = 5;
rps = 15;
kb = 10;

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
a = max(volfrac)/min(volfrac);

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

% round to nearest decimal
data(:,1) = round(data(:,1),1);

% find interval averages
r = intervals(boxstrain,sidex);
interval_average(data(:,1),data(:,4),r);

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



