clc;
close all;
clear;

%% Parameters
% Number of bins
File = fopen('Binning_parameters.txt','r');
data = fscanf(File,'%d',[3 1])';
fclose(File);

% Simulation parameters
nfib = 12800;
nseg = 5;
Lx = 1293;
rps = 15;
volfracAvg = nfib*2*pi*nseg*rps/(Lx*Lx*Lx); 

%% Array Allocation
nbin = data(1);
Bin = zeros(nbin,nbin,3);
BinHomo = zeros(nbin,nbin,3);
Box = zeros(nbin,nbin,nbin);
BoxHomo = zeros(nbin,nbin,nbin);

%% Read Files
% Flocculated xy, xz, yz
File = fopen('Binning_result_xy.txt','r');
data = fscanf(File,'%f',[nbin,nbin])';
fclose(File);
Bin(:,:,1) = data;

File = fopen('Binning_result_xz.txt','r');
data = fscanf(File,'%f',[nbin,nbin])';
fclose(File);
Bin(:,:,2) = data;

File = fopen('Binning_result_yz.txt','r');
data = fscanf(File,'%f',[nbin,nbin])';
fclose(File);
Bin(:,:,3) = data;

% Homogeneous xy, xz, yz
File = fopen('Binning_result_homo_xy.txt','r');
data = fscanf(File,'%f',[nbin,nbin])';
fclose(File);
BinHomo(:,:,1) = data;

File = fopen('Binning_result_homo_xz.txt','r');
data = fscanf(File,'%f',[nbin,nbin])';
fclose(File);
BinHomo(:,:,2) = data;

File = fopen('Binning_result_homo_yz.txt','r');
data = fscanf(File,'%f',[nbin,nbin])';
fclose(File);
BinHomo(:,:,3) = data;

% Flocculated Box
File = fopen('Binning_3d.txt','r');
data = fscanf(File,'%f',[nbin nbin*nbin])';
fclose(File);

% Reshape data
for i=1:nbin
    Box(:,:,i) = data((i-1)*nbin+1:i*nbin,:);
end

% Homogeneous Box
File = fopen('Binning_homo_3d.txt','r');
data = fscanf(File,'%f',[nbin nbin*nbin])';
fclose(File);

% Reshape data
for i=1:nbin
    BoxHomo(:,:,i) = data((i-1)*nbin+1:i*nbin,:);
end

%% Processing Test

%% Define Data To Be Analyzed
Name = {'\it{}xy\rm\bf{}','\it{}xz\rm\bf{}','\it{}yz\rm\bf{}'};
xtag = {'\it{}x','\it{}x','\it{}y'}; 
ytag = {'\it{}y','\it{}z','\it{}z'}; 

%% Mass Distribution
figure('color','white')
for i=1:3
    subplot(2,3,i)
    hold on
    contourf(Bin(:,:,i),'linestyle','none')
    title(['Segment C.M. (\mu = 20, ',Name{i}, ' view)'])
    xlabel(xtag{i}); 
    ylabel(ytag{i}); 
    colorbar
    box on
    set(gca,'fontsize',20,'fontName','Times New Roman','linewidth',2)
end

%% Intensity of Segregation (Based on Boxes)
Vcell = (Lx/nbin)^3;
Vseg = 2*pi*rps;
volfrac = Box*Vseg/Vcell;
volfracHomo = BoxHomo*Vseg/Vcell;
maxPhi = max(max(max(volfrac))) % check if volume fraction exceeds 1
if maxPhi > 1
    disp('Maximum \phi exceeds 1')
    pause
end
var = (volfrac-volfracAvg).^2;
varHomo = (volfracHomo-volfracAvg).^2;
I = mean(mean(mean(var)))/(volfracAvg*(1-volfracAvg));
IHomo = mean(mean(mean(varHomo)))/(volfracAvg*(1-volfracAvg));

%% Autocorrelation
avg = nfib*nseg/nbin^3;
auto = zeros(nbin/2,3);
autoHomo = zeros(nbin/2,3);
norm2Box = mean(mean(mean((Box-avg).^2)));
norm2BoxHomo = mean(mean(mean((BoxHomo-avg).^2)));
for offset=0:nbin/2-1
    sumx = 0;
    sumxHomo = 0;
    sumy = 0;
    sumyHomo = 0;
    sumz = 0;
    sumzHomo = 0;      
     for kk=0:nbin-1
        ind = kk+1;
        pairInd = kk+offset-floor((kk+offset)/nbin)*nbin+1;
        for j=1:nbin
            for i=1:nbin
                sumy = sumy + (Box(ind,i,j)-avg)*(Box(pairInd,i,j)-avg);
                sumyHomo = sumyHomo + (BoxHomo(ind,i,j)-avg)*(BoxHomo(pairInd,i,j)-avg);
                sumx = sumx + (Box(i,ind,j)-avg)*(Box(i,pairInd,j)-avg);
                sumxHomo = sumxHomo + (BoxHomo(i,ind,j)-avg)*(BoxHomo(i,pairInd,j)-avg);
                sumz = sumz + (Box(i,j,ind)-avg)*(Box(i,j,pairInd)-avg);
                sumzHomo = sumzHomo + (BoxHomo(i,j,ind)-avg)*(BoxHomo(i,j,pairInd)-avg);
            end
        end
    end
    auto(offset+1,1) = sumx;
    autoHomo(offset+1,1) = sumxHomo;
    auto(offset+1,2) = sumy;
    autoHomo(offset+1,2) = sumyHomo;
    auto(offset+1,3) = sumz;
    autoHomo(offset+1,3) = sumzHomo;
end
auto = auto/(nbin*nbin)/norm2Box;
autoHomo = autoHomo/(nbin*nbin)/norm2BoxHomo;
Direction = {'\it{}x\rm\bf{}','\it{}y\rm\bf{}','\it{}z\rm\bf{}'};
R = (0:nbin/2-1)/nbin*2;
auto = auto/nbin;
autoHomo = autoHomo/nbin;
hold on
for i=1:3
    subplot(2,3,i+3);
    hold on;
    box on;
    plot(R,auto(:,i),'color',rgb('MediumBlue'),'linewidth',2)
    plot(R,autoHomo(:,i),'--','color',rgb('Crimson'),'linewidth',2)
    title([Direction{i}, '-direction'])
    set(gca,'fontsize',20,'fontName','Times New Roman','linewidth',2)
    legend('\mu = 20','Frictionless')
	xlim([0 1])
    ylim([-0.2 1.01])
    xlabel('2\it{}r / L_{box}')
    ylabel('\it{}R(r)')
    set(gca,'YMinorTick','on','XMinorTick','on')
end

%% Length Scale
cross = (nbin/2-1)*ones(3,1);
crossHomo = (nbin/2-1)*ones(3,1);
l = zeros(3,1);
lHomo = zeros(3,1);
v = zeros(3,1);
vHomo = zeros(3,1);
autoV = zeros(nbin/2,3);
autoVHomo = zeros(nbin/2,3);
for j=1:3
    for i=1:nbin/2-1
        if auto(i,j) < 0
            cross(j) = i-1;
            break;
        end
    end
    for i=1:nbin/2-1
        if autoHomo(i,j) < 0
            crossHomo(j) = i-1;
            break;
        end
    end
end
R2 = R.^2;
R2 = R2';
for j=1:3
    autoV(:,j) = R2.*auto(:,j);
    autoVHomo(:,j) = R2.*autoHomo(:,j);
    l(j) = trapz(R(1:cross(j)),auto(1:cross(j),j));
    lHomo(j) = trapz(R(1:crossHomo(j)),autoHomo(1:crossHomo(j),j));
    v(j) = 2*pi*trapz(R(1:cross(j)),autoV(1:cross(j),j));
    vHomo(j) = 2*pi*trapz(R(1:crossHomo(j)),autoVHomo(1:crossHomo(j),j));
end

cross
crossHomo
S = l'
SHomo = lHomo'
V = v'
VHomo = vHomo'
