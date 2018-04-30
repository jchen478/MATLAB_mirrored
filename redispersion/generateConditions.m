clc;
clear;
close all;

%% input from file
% initial box
sidex = 600;
sidey = 600;
sidez = 600; 
% max conc = a * init conc
a = 2;
% volume change by rate*100% in concentrating
rate = 0.0005; 
% write info
box_write=1000;
dt=0.0001;

%% calculate strain in each region
N = ceil(-log(a)/log(1-rate)*dt*box_write); % strain in concentrating or expanding
region1 = 50;           % low conc steady shear
region2 = region1+N;  % concentration
region3 = region2+400;  % high conc steady shear
region4 = region3+N;  % expand
region5 = region4+200;  % low conc steady shear
gammaTot = region5; 


%% allocate array based on needed strain to complete cycle
gamma = (0:dt*box_write:gammaTot)';
Lx = zeros(length(gamma),1);
Ly = zeros(length(gamma),1);
Lz = zeros(length(gamma),1);

n = 0; 
%% fill in box array
for s=1:length(gamma)+1
    % region 1
    if (s <= region1/(dt*box_write)+1)
        Lx(s) = sidex;
        Ly(s) = sidey;
        Lz(s) = sidez;
        n = 1; 
        continue;
    end
    % region 2
    if (s <= region2/(dt*box_write)+1)
        Lx(s) = sidex*(1-rate)^(n/3);
        Ly(s) = sidey*(1-rate)^(n/3);
        Lz(s) = sidez*(1-rate)^(n/3);
        n = n+1; 
        continue;
    end
    % region 3
    if (s <= region3/(dt*box_write)+1)
        Lx(s) = Lx(s-1);
        Ly(s) = Ly(s-1);
        Lz(s) = Lz(s-1);
        ind = s; 
        n = 1; 
        continue;
    end
    % region 4
    if (s <= region4/(dt*box_write)+1)
        Lx(s) = Lx(ind)*(1-rate)^(-n/3);
        Ly(s) = Ly(ind)*(1-rate)^(-n/3);
        Lz(s) = Lz(ind)*(1-rate)^(-n/3);
        n = n+1; 
        continue;
    end
    % region 5
    if (s <= region5/(dt*box_write)+1)
        Lx(s) = sidex;
        Ly(s) = sidey;
        Lz(s) = sidez;
        continue;
    end
end
V = Lx.*Ly.*Lz;

figure()
hold on
plot(gamma, Lx)
plot(gamma, Ly)
plot(gamma, Lz)
legend('$L_x$','$L_y$','$L_z$')
xlabel('$\gamma$')
ylabel('Box length')

figure()
hold on
plot(gamma, V/(sidex*sidey*sidez))
xlabel('$\gamma$')
ylabel('$V/V_{original}$')
% 
% gammaTot