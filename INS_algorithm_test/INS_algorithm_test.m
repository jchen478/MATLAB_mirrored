%%{
clc
clear
close all
%}

%%{
% define common parameters / variables
side = 300;
nbin = 14;
rps = 15;
nsegTot = 800;
%}

name = 'rAscii.in';
File = fopen(name,'r');
for trials=1:5
    
    %%{
    % read files
    data = fscanf(File,'%f',[5 800])';
    % delete index and assign arrays
    data(:,1:2) = [];
    rx = data(:,1);
    ry = data(:,2);
    rz = data(:,3);
    %}
    
    %%{
    % shift coordinates so that all is in box
    shift = round(rx/side);
    rx = rx - side*shift;
    rx = rx + side/2;
    shift = round(ry/side);
    ry = ry - side*shift;
    ry = ry + side/2;
    shift = round(rz/side);
    rz = rz - side*shift;
    rz = rz + side/2;
    % Visualize fiber distribution
        figure()
        scatter3(rx,ry,rz,'filled','markerFaceColor',rgb('MediumBlue'))
        xlabel('$x$')
        ylabel('$y$')
        zlabel('$z$')
    %}
    
    I = intensityCalc( rx, ry, rz, side, nbin, nsegTot, rps);
end
fclose(File);

close all;

%%{
% Test cases
factorArr = 0.3:0.05:1;
% factorArr = 1;
nFactor = length(factorArr);
nbinArr = [17 18 20 22 24 26 30 32 34 40];
nNbin = length(nbinArr);
nbinLegendArr = cell(nNbin,1); 
IArr = zeros(nFactor,nNbin);
dbinArr = side ./nbinArr; 

for i=1:nFactor
    
    factor = factorArr(i);
    rx = factor*side*rand(nsegTot,1);
    ry = factor*side*rand(nsegTot,1);
    rz = factor*side*rand(nsegTot,1);
    
    for j=1:nNbin
        
        nbin = nbinArr(j);
        IArr(i,j) = intensityCalc( rx, ry, rz, side, nbin, nsegTot, rps);
        
        %     figure(i)
        %     scatter3(rx,ry,rz,'filled','markerFaceColor',rgb('MediumBlue'))
        %     xlabel('$x$')
        %     ylabel('$y$')
        %     zlabel('$z$')
        %     xlim([0 side])
        %     ylim([0 side])
        %     zlim([0 side])
        %     title(['$I = $ ', num2str(IArr(i),'%.3f')])
        
        %}
    end
end

figure(nFactor + 1)
hold on
for j=1:nNbin
    
%     plot(factorArr,IArr(:,j),'-.o','markerFaceColor',rgb('MediumBlue'),...
%         'markerEdgeColor',rgb('MediumBlue'),'linewidth',2.5,'color',rgb('MediumBlue'))
    plot(factorArr,IArr(:,j),'-.o','linewidth',2.5)
    xlabel('Fraction occupied $L_{box}$')
    ylabel('$I$')
    nbinLegendArr{j} = ['$N_{bin} =$ ',num2str(nbinArr(j)),' $d_{bin} =$ ', num2str(dbinArr(j),'%.1f')]; 
    
end
xlim([0 1])
legend(nbinLegendArr,'location','best')