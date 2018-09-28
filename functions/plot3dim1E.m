function [] = plot3dim1E(data1, data1E, x, y, z)
%PLOT3DIM Summary of this function goes here
%   Given 3d data and x-axis and legend of interest

%% permute data according to dimension
pdata1 = permute(data1.value,[x.dim y.dim z.dim]);
pdata2 = permute(data1E,[x.dim y.dim z.dim]);

%% create corresponding legend array
for k=1:length(z.value)
    figure('Units','Inches','Position',[4 6 3.0 2.5]);
%     figure('Units','Inches','Position',[1 1 3.5 3.0]);
    hold on
    set(gcf,'name',z.legend{k})
    for j=1:length(y.value)
        errorbar(x.value,pdata1(:,j,k),pdata2(:,j,k),'-.o','MarkerSize',6, 'linewidth',1.8)
    end
    box on
    hold on
    title(z.legend{k})
    xlabel(x.name)
    ylabel(data1.name)
    legend(y.legend,'location','best')  
end
end

