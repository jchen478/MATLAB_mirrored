function elasticFun(filename,Seff)

%% Elastic storage
File = fopen(filename,'r');
elastic = fscanf(File,'%f',[4 Inf])';
fclose(File);



%% Figure formatting parameters
markersize = 10;
tickx = 1.5;
fontsize = 24;
linewidth = 1.5;


%% Plotting
% figure()
hold on
scatter(elastic(:,1)*Seff*(410/pi),elastic(:,4)/elastic(1,4),...
        markersize,'MarkerEdgeColor',rgb('MediumBlue'),...
        'MarkerFaceColor',rgb('MediumBlue'))
ylabel('\it{}E_{el}')



end
