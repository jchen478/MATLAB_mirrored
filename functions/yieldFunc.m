function [yieldstress] = yieldFunc( filename, ID, yieldstress, nconc )


    %% read file
    File = fopen(filename,'r');
    data = fscanf(File,'%f',[11 Inf])';
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
    color = {rgb('Maroon'),rgb('DeepPink'),rgb('DarkOrange'),rgb('DarkGreen'),rgb('MediumBlue'),rgb('DarkMagenta')};

    regressX = linspace(0,max(over_Seff),400);
    j = 1;
    figure(ID)
    hold on;
    for i=1:3:length(kb)
        
        
        s(j) = errorbar(over_Seff(i:i+2),sigxzStar(i:i+2),...
            se_sigxzStar(i:i+2)/2,'o','MarkerSize',10,'MarkerEdgeColor',color{j},'Linewidth',2.5,'color',color{j});
%         s(j) = scatter(over_Seff(i:i+2),sigxzStar(i:i+2),100,'MarkerEdgeColor',color{j},'Linewidth',2.0);
        over_Seff(i:i+2);
        sigxzStar(i:i+2);
        fit = polyfit(over_Seff(i:i+2),sigxzStar(i:i+2),1);
        sigxz_regress = fit(1)*regressX+fit(2);
        h(j) = plot(regressX,sigxz_regress,'color',color{j},'linewidth',2);
        sig0(j) = fit(2);
        sig0_phi(j) = volfrac(i); 
        S{2*j} = ['\sigma_{xz} = ', num2str(fit(1),3), ' / S_{eff} + ', num2str(fit(2),'%.2f')];
        j = j+1;

        
    end
    title(filename)
            
    yieldstress((ID-1)*nconc+1:nconc*ID,1) = nfib(1:nconc);
    yieldstress((ID-1)*nconc+1:nconc*ID,2) = sig0_phi;
    yieldstress((ID-1)*nconc+1:nconc*ID,3) = sig0;
    
end

