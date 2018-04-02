function [ ] = evolution_stress_changing_mu(name,color,nfib,nseg,lbox,rps,kb,a,EY,Imom,eta0)
%   plot stress 

%% open file
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

%% calculate relevant parameters
nfibArray = nfib*ones(nStep,1);
nsegArr = nseg*ones(nStep,1);
side = lbox*ones(nStep,1);
[Seff,rp,nL3,L,gamma] = pCalc(nfibArray, nsegArr, rps, kb, side, a, EY, Imom, eta0);

% dimensionalize stresses
[sigmap,sigma, std_sigmap, std_sigma, N1, N2] = ...
    stressDim(data(:,4), zeros(nStep,1), zeros(nStep,1), ...
    zeros(nStep,1), nseg, rps, eta0, gamma, nL3);

% nondimensionalize using natural variables 
sigmap_nondim = sigmap.*L.^4/EY/Imom;
strain = data(:,1);

plot(strain,sigmap_nondim,'-o','color',color,...
    'MarkerFaceColor',color,...
    'MarkerEdgeColor',color,...
    'MarkerSize',2)
end

