close all
clc
clear

%% Figure formatting parameters
nfig = 12;
tickx = 1.5;
fontsize = 24;
linewidth = 3;
File = fopen('test.txt','r');
data = fscanf(File,'%f',[4 Inf])';
fclose(File);

%% Parameters
L = 150;

dx = 30;
dy = 30;
dz = 30;

% sidex = 600;
% sidey = 600;
% sidez = 600;

sidex = 450;
sidey = 450;
sidez = 1066;


rx = data(:,2)-sidex;
ry = data(:,3)-sidey;
rz = data(:,4)-sidez;

N = length(rx);

samplex = 10;
sampley = 10;
samplez = 22;
maxGap = 2;

%% Visualize cluster
figure(1)
subplot(1,2,1)
hold on
box on
scatter3(rx,ry,rz)
xlim([-sidex/2 sidex/2])
ylim([-sidey/2 sidey/2])
zlim([-sidez/2 sidez/2])
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',tickx*get(gca,'ticklength'))
set(gcf, 'color','white')
title('Original')

%% Shift segments outside of box into box
for i=1:N
    if (rx(i) < -sidex/2) || (rx(i) > sidex/2)
        rx(i) = rx(i) - sign(rx(i))* sidex;
    end
    if (ry(i) < -sidey/2) || (ry(i) > sidey/2)
        ry(i) = ry(i) - sign(ry(i))* sidey;
    end
    if (rz(i) < -sidez/2) || (rz(i) > sidez/2)
        rz(i) = rz(i) - sign(rz(i))* sidez;
    end
end

%% Shift segments so that all coordinates are positive
rx = rx+sidex/2;
ry = ry+sidey/2;
rz = rz+sidez/2;

%% Visualize shifted cluster
figure(1)
subplot(1,2,2)
hold on
box on
scatter3(rx,ry,rz,'filled','markerFaceColor',rgb('Crimson'))
xlabel('x')
ylabel('y')
zlabel('z')
xlim([0 sidex])
ylim([0 sidey])
zlim([0 sidez])
set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',tickx*get(gca,'ticklength'))
set(gcf, 'color','white')
title('Shifted')

%% Divide domain into bins
nx = ceil(sidex/dx);
ny = ceil(sidey/dy);
nz = ceil(sidez/dz);
bin = zeros(ny,nx,nz);

dx = sidex/nx;
dy = sidey/ny;
dz = sidez/nz;

%% Place fiber segments into bin
for i=1:N
    xloc = ceil(rx(i)/dx);
    yloc = ceil(ry(i)/dy);
    zloc = ceil(rz(i)/dz);
    bin(yloc,xloc,zloc) =  bin(yloc,xloc,zloc) + 1;
end


%% Visualize binning results
k = 1;
A = zeros(ny,nx);
for i=1:nx
    for j=1:ny
        A(j,i) = bin(j,i,k);
    end
end

figure(2)
subplot(1,2,1)
imagesc(A)
box on
xlabel('Bin_x')
ylabel('Bin_y')
title('Before filling x')
set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',tickx*get(gca,'ticklength'))
set(gcf, 'color','white')


%% Fill in blanks in x direction
binx = zeros(ny,nx,nz);
for j=1:ny
    for k = 1:nz
        if(max(bin(j,:,k)) > 0)
            for i = nx:-1:1
                num = bin(j,i,k);
                if (num ~= 0)
                    oldBin = i;
                    break;
                end
            end
            newBin = 1;
            for i = 1:nx
                num = bin(j,i,k);
                if (num ~= 0)
                    binx(j,i,k) = 1;
                    newBin = i;
                    if(oldBin > newBin)
                        if(nx-oldBin+newBin <= maxGap)
                            for ii=oldBin:nx
                                binx(j,ii,k) = 1;
                            end
                            for ii=1:newBin
                                binx(j,ii,k) = 1;
                            end
                        end
                    else
                        if(newBin-oldBin <= maxGap)
                            for ii=oldBin:newBin
                                binx(j,ii,k) = 1;
                            end
                        end
                    end
                    oldBin = newBin;
                end
            end
        end
    end
end

%% Visualize binning results
k = 1;
A = zeros(ny,nx);
for i=1:nx
    for j=1:ny
        A(j,i) = binx(j,i,k);
    end
end

figure(2)
subplot(1,2,2)
imagesc(A)
box on
xlabel('Bin_x')
ylabel('Bin_y')
title('After filling x')
set(gca,'fontsize',fontsize,'linewidth',linewidth,'fontname','Times New Roman')
set(gca,'YMinorTick','on','XMinorTick','on')
set(gca,'ticklength',tickx*get(gca,'ticklength'))
set(gcf, 'color','white')

%% Calculate mean x-length
xdim = 0;
ind = 1;
xdimCutoff = 1;
for j=1:ny
    for k=1:nz
        if (sum(binx(j,:,k)) > xdimCutoff)
            xdim(ind) = sum(binx(j,:,k));
            ind = ind+1;
        end
    end
end

xdim = xdim';

xdim = max(xdim)*dx

%% Fill in blanks in y direction
biny = zeros(ny,nx,nz);
for i=1:nx
    for k = 1:nz
        if(max(bin(:,i,k)) > 0)
            for j = ny:-1:1
                num = bin(j,i,k);
                if (num ~= 0)
                    oldBin = j;
                    break;
                end
            end
            newBin = 1;
            for j = 1:ny
                num = bin(j,i,k);
                if (num ~= 0)
                    biny(j,i,k) = 1;
                    newBin = j;
                    if(oldBin > newBin)
                        if(ny-oldBin+newBin <= maxGap)
                            for jj=oldBin:ny
                                biny(jj,i,k) = 1;
                            end
                            for jj=1:newBin
                                biny(jj,i,k) = 1;
                            end
                        end
                    else
                        if(newBin-oldBin <= maxGap)
                            for jj=oldBin:newBin
                                biny(jj,i,k) = 1;
                            end
                        end
                    end
                    oldBin = newBin;
                end
            end
        end
    end
end

%% Calculate mean y-length
ydim = 0;
ind = 1;
ydimCutoff = 1;
for i=1:nx
    for k=1:nz
        if (sum(biny(:,i,k)) > ydimCutoff)
            ydim(ind) = sum(biny(:,i,k));
            ind = ind+1;
        end
    end
end

ydim = ydim';

ydim = max(ydim)*dy


%% Fill in blanks in z direction
binz = zeros(ny,nx,nz);
for i=1:nx
    for j = 1:ny
        access = false;
        if(max(bin(j,i,:)) > 0)
            for k = nz:-1:1
                num = bin(j,i,k);
                if (num ~= 0)
                    oldBin = k;
                    access = true;
                    break;
                end
            end
            newBin = 1;
            for k = 1:nz
                num = bin(j,i,k);
                if (num ~= 0)
                    binz(j,i,k) = 1;
                    newBin = k;
                    if(access == true)
                        if(nz-oldBin+newBin <= maxGap)
                            for kk=oldBin:nz
                                binz(j,i,kk) = 1;
                            end
                            for kk=1:newBin
                                binz(j,i,kk) = 1;
                            end
                        end
                        access = false;
                    else
                        if(newBin-oldBin <= maxGap)
                            for kk=oldBin:newBin
                                binz(j,i,kk) = 1;
                            end
                        end
                    end
                    oldBin = newBin;
                end
            end
        end
    end
end

%% Calculate mean z-length
zdim = 0;
ind = 1;
zdimCutoff = 1;
for i=1:nx
    for j=1:ny
        if (sum(binz(j,i,:)) > zdimCutoff)
            zdim(ind) = sum(binz(j,i,:));
            ind = ind+1;
        end
    end
end

zdim = zdim';

zdim = max(zdim)*dz
