function [] = plot2dim1EB(data1, data1E, x, y)
%PLOT3DIM Summary of this function goes here
%   Given 3d data and x-axis and legend of interest

%% permute data according to dimension
pdata1 = permute(data1.value,[x.dim y.dim]);
pdata2 = permute(data1E,[x.dim y.dim]);

%% create corresponding legend array
figure('Units','Inches','Position',[1 1 3.0 2.5]);
hold on
for j=1:length(y.value)
    errorbar(x.value,pdata1(:,j),pdata2(:,j),'-.o','MarkerSize',6, 'linewidth',1.8)
end
box on
hold on
xlabel(x.name)
ylabel(data1.name)
legend(y.legend,'location','best')

end


