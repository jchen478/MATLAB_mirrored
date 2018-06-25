function [] = plotPhase(data1, data1Name, x, y, z, criteria)
%PLOT3DIM Summary of this function goes here
%   Given 3d data and x-axis and legend of interest

%% permute data according to dimension
pdata1 = permute(data1,[x.dim y.dim z.dim]);

%% create corresponding legend array
for k=1:length(z.value)
    figure('Units','Inches','Position',[1 1 3.2 3.0]);
    set(gcf,'name',z.legend{k})
    for j=1:length(y.value)
        
        % convert data to logicals based on criteria
        dataRaw = pdata1(:,j,k);
        dataBelow = dataRaw < criteria;
        dataAbove = dataRaw >= criteria;
        
        % convert logicals to double and then zero elements to NaN
        dataBelow = double(dataBelow);
        dataBelow(dataBelow == 0) = nan;
        dataAbove = double(dataAbove);
        dataAbove(dataAbove == 0) = nan;
        
        box on
        hold on
        scatter(x.value,y.value(j)*dataBelow,30,'o',...
            'MarkerFaceColor',rgb('Black'),...
            'MarkerEdgeColor',rgb('Black'),...
            'LineWidth',2.0)
        scatter(x.value,y.value(j)*dataAbove,70,'x',...
            'MarkerEdgeColor',rgb('Black'),...
            'LineWidth',2.0)
        if j == 1
            legendArr = {[data1Name,' $< $ ',num2str(criteria)],...
                [data1Name,' $\geq $ ',num2str(criteria)]};
            legend(legendArr,'location','best')
        end
    end
%     title(z.legend{k})
    xlabel(x.name)
    ylabel(y.name)
end
end

