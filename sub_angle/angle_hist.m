clc
clear
% close all

%% histogram definition
binWidth = 0.05;
nBin = 32;
angLim = nBin*binWidth;

%% bending histogram
filename = 'dBend.txt';
File = fopen(filename,'r');
dBend = fscanf(File,'%f',[33 Inf])';
fclose(File);
strain = dBend(:,1);
strain = round(strain,1); 
dBend(:,1) = [];
nStrain = length(strain);

%% twisting histogram
filename = 'dTwist.txt';
File = fopen(filename,'r');
dTwist = fscanf(File,'%f',[33 Inf])';
fclose(File);
dTwist(:,1) = [];

%% angle statistics
filename = 'Angle.txt';
File = fopen(filename,'r');
AngleStat = fscanf(File,'%f',[7 Inf])';
fclose(File);

%% plot
figure()
hold on
maxBend = max(max(dBend(790*2:1090*2,:)));
for i=790*2:1090*2
    plot(binWidth/2:binWidth:angLim-binWidth/2,dBend(i,:),'-o')
    ylim([0 maxBend]);
%     M(i) = getframe(gcf);
end

%%
% figure()
% axes('Position',[0 0 1 1])
% movie(M,1)