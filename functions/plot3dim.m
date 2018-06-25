function [] = plot3dim(data1, data2, data1Name, data2Name, x, y, z)
%PLOT3DIM Summary of this function goes here
%   Given 3d data and x-axis and legend of interest

%% permute data according to dimension
pdata1 = permute(data1,[x.dim y.dim z.dim]);
pdata2 = permute(data2,[x.dim y.dim z.dim]);

%% create corresponding legend array
for k=1:length(z.value)
    figure('Units','Inches','Position',[1 1 4.0 6.0]);
    set(gcf,'name',z.legend{k})
    for j=1:length(y.value)
        
        subplot(2,1,1)
        hold on
        plot(x.value,pdata1(:,j,k),'-o')

        subplot(2,1,2)
        hold on
        plot(x.value,pdata2(:,j,k),'-o')
        
    end
    
    subplot(2,1,1)
    box on
    hold on
    title(z.legend{k})
    xlabel(x.name)
    ylabel(data1Name)
    legend(y.legend,'location','best')
    
    subplot(2,1,2)
    box on
    hold on
%     title('\bf{Number of contacts}')
    xlabel(x.name)
    ylabel(data2Name)
    legend(y.legend,'location','best')
    
end
end

