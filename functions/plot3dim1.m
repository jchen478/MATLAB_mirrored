function [] = plot3dim1(data1, x, y, z)
%PLOT3DIM Summary of this function goes here
%   Given 3d data and x-axis and legend of interest

%% permute data according to dimension
pdata1 = permute(data1.value,[x.dim y.dim z.dim]);

%% create corresponding legend array
for k=1:length(z.value)
    figure('Units','Inches','Position',[1 1 3.0 2.5]);
    hold on
    set(gcf,'name',z.legend{k})
    for j=1:length(y.value)
        plot(x.value,pdata1(:,j,k),'-o')
    end
    box on
    hold on
    title(z.legend{k})
    xlabel(x.name)
    ylabel(data1.name)
    legend(y.legend,'location','best')  
end
end

