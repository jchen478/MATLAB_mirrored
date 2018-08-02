% clc;
% clear;
% close all;

%% input from file
% initial box
sidex = 600;
sidey = 600;
sidez = 600;

% max conc = a * init conc
a = 5;
% volume change by rate*100% in concentrating
rate = 0.0005;
% write info
box_write=1000;
dt=0.0001;

% specify expansion (-) or compression (+)
pow3 = 0.6;
pow1 = 0.2;
pow2 = 0.2;

%% Check binning
dxConst = 18.0;
dyConst = 18.0;
dzConst = 18.0;

nxbinMax = floor(sidex / dxConst);
nybinMax = floor(sidey / dyConst);
nzbinMax = floor(sidez / dzConst);

if (mod(nxbinMax,2) ~= 0)
    nxbinMax = nxbinMax - 1;
end
if (mod(nybinMax,2) ~= 0)
    nybinMax = nybinMax - 1;
end
if (mod(nzbinMax,2) ~= 0)
    nzbinMax = nzbinMax - 1;
end

%% calculate strain in each region
N = ceil(-log(a)/log(1-rate)*dt*box_write); % strain in concentrating or expanding
region1 = 100;           % low conc steady shear
region2 = region1+N;  % concentration
region3 = region2+300;  % high conc steady shear
region4 = region3+N;  % expand
region5 = region4+300;  % low conc steady shear
gammaTot = region5;


%% allocate array based on needed strain to complete cycle
gamma = (0:dt*box_write:gammaTot)';
Lx = zeros(length(gamma),1);
Ly = zeros(length(gamma),1);
Lz = zeros(length(gamma),1);

binx = zeros(length(gamma),1);
biny = zeros(length(gamma),1);
binz = zeros(length(gamma),1);

n = 0;
%% fill in box array
for s=1:length(gamma)+1
    % region 1
    if (s <= region1/(dt*box_write)+1)
        Lx(s) = sidex;
        Ly(s) = sidey;
        Lz(s) = sidez;
        n = 1;
        binx(s) = floor(Lx(s) / dxConst);
        biny(s) = floor(Ly(s) / dyConst);
        binz(s) = floor(Lz(s) / dzConst);
        if (mod(binx(s),2) ~= 0)
            binx(s) = binx(s) - 1;
        end
        if (mod(biny(s),2) ~= 0)
            biny(s) = biny(s) - 1;
        end
        if (mod(binz(s),2) ~= 0)
            binz(s) = binz(s) - 1;
        end
        continue;
    end
    % region 2
    if (s <= region2/(dt*box_write)+1)
        Lx(s) = sidex*(1-rate)^(n*pow1);
        Ly(s) = sidey*(1-rate)^(n*pow2);
        Lz(s) = sidez*(1-rate)^(n*pow3);
        n = n+1;
        binx(s) = floor(Lx(s) / dxConst);
        biny(s) = floor(Ly(s) / dyConst);
        binz(s) = floor(Lz(s) / dzConst);
        if (mod(binx(s),2) ~= 0)
            binx(s) = binx(s) - 1;
        end
        if (mod(biny(s),2) ~= 0)
            biny(s) = biny(s) - 1;
        end
        if (mod(binz(s),2) ~= 0)
            binz(s) = binz(s) - 1;
        end
        continue;
    end
    % region 3
    if (s <= region3/(dt*box_write)+1)
        Lx(s) = Lx(s-1);
        Ly(s) = Ly(s-1);
        Lz(s) = Lz(s-1);
        ind = s;
        n = 1;
        binx(s) = floor(Lx(s) / dxConst);
        biny(s) = floor(Ly(s) / dyConst);
        binz(s) = floor(Lz(s) / dzConst);
        if (mod(binx(s),2) ~= 0)
            binx(s) = binx(s) - 1;
        end
        if (mod(biny(s),2) ~= 0)
            biny(s) = biny(s) - 1;
        end
        if (mod(binz(s),2) ~= 0)
            binz(s) = binz(s) - 1;
        end
        continue;
    end
    % region 4
    if (s <= region4/(dt*box_write)+1)
        Lx(s) = Lx(ind)*(1-rate)^(-n*pow1);
        Ly(s) = Ly(ind)*(1-rate)^(-n*pow2);
        Lz(s) = Lz(ind)*(1-rate)^(-n*pow3);
        n = n+1;
        binx(s) = floor(Lx(s) / dxConst);
        biny(s) = floor(Ly(s) / dyConst);
        binz(s) = floor(Lz(s) / dzConst);
        if (mod(binx(s),2) ~= 0)
            binx(s) = binx(s) - 1;
        end
        if (mod(biny(s),2) ~= 0)
            biny(s) = biny(s) - 1;
        end
        if (mod(binz(s),2) ~= 0)
            binz(s) = binz(s) - 1;
        end
        continue;
    end
    % region 5
    if (s <= region5/(dt*box_write)+1)
        Lx(s) = sidex;
        Ly(s) = sidey;
        Lz(s) = sidez;
        binx(s) = floor(Lx(s) / dxConst);
        biny(s) = floor(Ly(s) / dyConst);
        binz(s) = floor(Lz(s) / dzConst);
        if (mod(binx(s),2) ~= 0)
            binx(s) = binx(s) - 1;
        end
        if (mod(biny(s),2) ~= 0)
            biny(s) = biny(s) - 1;
        end
        if (mod(binz(s),2) ~= 0)
            binz(s) = binz(s) - 1;
        end
        continue;
    end
end
V = Lx.*Ly.*Lz;

% figure()
% hold on
% plot(gamma, Lx)
% plot(gamma, Ly)
% plot(gamma, Lz)
% legend('$L_x$','$L_y$','$L_z$')
% xlabel('$\gamma$')
% ylabel('Box length')
% 
% figure()
% hold on
% plot(gamma, V/(sidex*sidey*sidez))
% xlabel('$\gamma$')
% ylabel('$V/V_{original}$')
%
% gammaTot

% figure
% hold on
% plot(gamma, binx)
% plot(gamma, biny)
% plot(gamma, binz)
% xlabel('$\gamma$')
% ylabel('$N_{bin}$')
% legend('$N_x$','$N_y$','$N_z$')
% 
% figure
% hold on
% plot(gamma, binx.*biny.*binz/(nxbinMax*nybinMax*nzbinMax))
% xlabel('$\gamma$')
% ylabel('$N_{bin}^3 / N_{bin,max}^3$')

if max(binx.*biny.*binz/(nxbinMax*nybinMax*nzbinMax)) > 1
    display('Error')
else
    display('This condition works')
end

volfrac = 1280*150*pi./(Lx.*Ly.*Lz);

% figure('Units','Inches','Position',[1 1 3.5 2.2]);
hold on
plot(gamma, volfrac*100)
% plot(gamma, volfrac*100,'color',rgb('black'))
% legend('$r_{\phi}=4$')
xlim([0 max(gamma)])
xlabel('$\gamma$')
ylabel('$\phi\ \%$')

% figure('Units','Inches','Position',[1 1 5.0 2.5]);
% hold on
% plot(gamma, volfrac*100,'color',rgb('black'))
% % legend('$r_{\phi}=4$')
% xlim([0 max(gamma)])
% xlabel('$\gamma$')
% ylabel('$\phi\ \%$')
% 
