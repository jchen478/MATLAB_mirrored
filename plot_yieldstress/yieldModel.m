function [] = yieldModel( filename, ID, yieldstress, nfit, nstiff, color)


%% read file
File = fopen(filename,'r');
data = fscanf(File,'%f',[26 Inf])';
fclose(File);

%% assign variable
nfib = data(:,1);
nseg = data(:,2);
rps = data(:,3);
kb = data(:,4);
mu = data(:,5);
side = data(:,6);
tauxzstar = data(:,7);
N1star = data(:,8);
N2star = data(:,9);
nc = data(:,10);
se_tauxzstar = data(:,11);

%% calculate parameters
L = 2*rps.*nseg;
nL3 = nfib.*L.^3./side.^3;
Seff = kb*pi./nseg.^4;
over_Seff = 1./Seff;
volfrac = pi*L.*nfib./side.^3;

%% total suspension stress
sigxzStar = tauxzstar*pi.*nL3./(6*nseg.^3.*log(2*rps))./Seff + 1./Seff;
se_sigxzStar = se_tauxzstar*pi.*nL3./(6*nseg.^3.*log(2*rps))./Seff;

%% plotting

regressX = linspace(0,max(over_Seff),400);
hold on;
for i=1:nstiff:length(kb)
    
    %% Bingham fit
    %  Fit
    figure(1)
    hold on;
    errorbar(over_Seff(i:i+nstiff-1),sigxzStar(i:i+nstiff-1),...
        se_sigxzStar(i:i+nstiff-1)/2,'o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
    x = over_Seff(i:i+nstiff-1);
    y = sigxzStar(i:i+nstiff-1);
    fit = polyfit(x,y,1);
    sig0_Bingham = fit(2)
    fit(1)
    sigxz_regress = fit(1)*regressX+fit(2);
    plot(regressX,sigxz_regress,'color',color,'linewidth',2);
    title('Bingham fit')
    % Error
    yFit = fit(1)*x+fit(2);
    RSS = (y-yFit).^2;
    figure(2)
    hold on;
    plot(x,RSS,'-o','linewidth',2,'markersize',10,'color',color,'color',color);
    title('Bingham fit residual error squared')
    % Error sum
    RSS_Bingham = sum(RSS)
    
    %% Casson fit
    %  Fit
    figure(3)
    hold on;
    errorbar(over_Seff(i:i+nstiff-1),sigxzStar(i:i+nstiff-1),...
        se_sigxzStar(i:i+nstiff-1)/2,'o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
    x = sqrt(over_Seff(i:i+nstiff-1));
    y = sqrt(sigxzStar(i:i+nstiff-1));
    fit = polyfit(x,y,1);
    sigxz_regress = fit(1)*sqrt(regressX)+fit(2);
    sigxz_regress = sigxz_regress.^2;
    plot(regressX,sigxz_regress,'color',color,'linewidth',2);
    title('Casson fit')
    % Error
    yFit = fit(1)*x+fit(2);
    sig0_Casson = fit(2)^2
    fit(1)^2
    yFit = yFit.^2;
    y = y.^2;
    RSS = (y-yFit).^2;
    figure(4)
    hold on;
    plot(x.^2,RSS,'-o','linewidth',2,'markersize',10,'color',color,'color',color);
    title('Casson fit residual error squared')
    % Error sum
    RSS_Casson = sum(RSS)
    
    %% Herschel-Bulkley fit
    %  Fit
    figure(5)
    hold on;
    errorbar(over_Seff(i:i+nstiff-1),sigxzStar(i:i+nstiff-1),...
        se_sigxzStar(i:i+nstiff-1)/2,'o','MarkerSize',10,'MarkerEdgeColor',color,'Linewidth',2.5,'color',color);
    x = over_Seff(i:i+nstiff-1);
    y = sigxzStar(i:i+nstiff-1);
    if ID == 1
        sig0 = 3.173729193;
        k = 28.88249607;
        n = 0.807807105;
    elseif ID == 2
        sig0 = 0;
        k = 25.54492471;
        n = 0.622483434;
    else
        sig0 = 2.808516336;
        k = 2.006747949;
        n = 0.928995207;
    end
    
            sig0 = 2.970898498
        k = 3.258936603;
        n = 0.71823196;
    
    sigxz_regress = k*regressX.^n+sig0;
    plot(regressX,sigxz_regress,'color',color,'linewidth',2);
    title('Herschel-Bulkley fit')
    % Error
    yFit = k*x.^n+sig0;
    RSS = (y-yFit).^2;
    figure(6)
    hold on;
    plot(x,RSS,'-o','linewidth',2,'markersize',10,'color',color,'color',color);
    title('Herschel-Bulkley fit residual error squared')
    % Error sum
    RSS_HB = sum(RSS)
end

end

