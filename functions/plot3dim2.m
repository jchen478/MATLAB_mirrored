function [] = plot3dim2(data1, data2, x, y, z)
%PLOT3DIM Summary of this function goes here
%   Given 3d data and x-axis and legend of interest

%% permute data according to dimension
pdata1 = permute(data1.value,[x.dim y.dim z.dim]);
pdata2 = permute(data2.value,[x.dim y.dim z.dim]);

%% create corresponding legend array
for k=1:length(z.value)
    figure('Units','Inches','Position',[1 1 3.0 5.5]);
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
    ylabel(data1.name)
    legend(y.legend,'location','best')
    
    subplot(2,1,2)
    box on
    hold on
    xlabel(x.name)
    ylabel(data2.name)
    legend(y.legend,'location','best')
    
end
end

