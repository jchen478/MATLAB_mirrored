function yieldStress(stress,title_nfib,loc)

% nfib nL3 kb volfrac sigxz stdev sem contact 1 2

sigxz_p = stress(:,5);
kb = stress(:,3); 
volfrac = stress(:,4);
nL3 = stress(:,2); 
SEM = stress(:,7);

nseg = 5;
rps = 15;
EYI = 1.128e-10; % N m2
L = 3.45e-3; % m

Seff = pi*kb/nseg^4; 
over_Seff = 1./Seff;

sigxz_total = sigxz_p*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff + over_Seff;
SEM = SEM*pi.*nL3/(6*nseg^3*log(2*rps)).*over_Seff/2; 

j = 0;
S = cell(2*length(kb)/3,1);
Seta = cell(2*length(kb)/3,1);
sig0 = zeros(length(kb)/3,1); 
sig0_phi = zeros(length(kb)/3,1); 

for i=1:3:length(kb)
    
    S{2*j+1} = ['\phi = ', num2str(volfrac(i),'%.3f')];
    Seta{2*j+1} = ['\phi = ', num2str(volfrac(i),'%.3f')];
    j = j+1; 
    
end

color = {rgb('DarkRed'),rgb('DarkOrange'),rgb('MediumBlue'),rgb('DarkOliveGreen')};

figure(1)
% subplot(2,2,loc)
box on;
hold on;

figure(3)
box on;
hold on;

regressX = linspace(0,max(over_Seff),400);
j = 1; 

for i=1:3:length(kb)
    
  
    figure(1)
    hold on;
    scatter(over_Seff(i:i+2),sigxz_total(i:i+2),40,'MarkerEdgeColor',color{j},...
              'MarkerFaceColor',color{j}); 
    fit = polyfit(over_Seff(i:i+2),sigxz_total(i:i+2),1);
    sigxz_regress = fit(1)*regressX+fit(2); 
    plot(regressX,sigxz_regress,'color',color{j},'linewidth',1.3);
    sig0(j) = fit(2); 
    sig0_phi(j) = volfrac(i); 
    S{2*j} = ['\sigma_{xz} = ', num2str(fit(1),3), ' / S_{eff} + ', num2str(fit(2),3)];
    
    figure(3)
    hold on;
    scatter(over_Seff(i:i+2),sigxz_total(i:i+2)./over_Seff(i:i+2),40,'MarkerEdgeColor',color{j},...
              'MarkerFaceColor',color{j}); 
    fit = polyfit(over_Seff(i:i+2),sigxz_total(i:i+2)./over_Seff(i:i+2),1);
    sigxz_regress = fit(1)*regressX+fit(2); 
    plot(regressX,sigxz_regress,'color',color{j},'linewidth',1.3);
    sig0(j) = fit(2); 
    sig0_phi(j) = volfrac(i); 
    Seta{2*j} = ['\eta* = ', num2str(fit(1),3), ' / S_{eff} + ', num2str(fit(2),3)];
    j = j+1; 
    
    sigxz_total(i:i+2)./over_Seff(i:i+2)
    
end

figure(1)
hold on
xlabel('1/S_{eff}')
ylabel('\sigma_{xz} L^4/E_YI')
set(gca,'fontsize',14,'fontName','Times New Roman')
set(gcf,'color','white')
legend(S{:})

figure(3)
hold on
xlabel('1/S_{eff}')
ylabel('\eta / \eta_0')
set(gca,'fontsize',36,'fontName','Times New Roman')
set(gcf,'color','white')
% legend(Seta{:})
legend({'\it{}\Phi\rm{} = 0.013',...
' \it{}\eta*\rm{} = (-0.353\pm0.054) / \it{}S_{eff}\rm{} + (22.35\pm1.43)',...
'\it{}\Phi\rm{} = 0.007',...
' \it{}\eta*\rm{} = (-0.025\pm0.001) / \it{}S_{eff}\rm{} + (2.50\pm0.05)',...
'\it{}\Phi\rm{} = 0.004',...
' \it{}\eta*\rm{} = (-0.00165\pm0.00002) / \it{}S_{eff}\rm{} + (1.2127\pm0.0005)',...
'\it{}\Phi\rm{} = 0.003',...
' \it{}\eta*\rm{} = (-0.00005\pm0.00015) / \it{}S_{eff}\rm{} + (1.070\pm0.004)'},...
'Location','northeast');



if(sig0(4) < 0)

figure(2)
subplot(2,2,loc)
hold on;
scatter(sig0_phi(1:3),sig0(1:3),40,'MarkerEdgeColor',rgb('MediumBlue'),...
              'MarkerFaceColor',rgb('MediumBlue'));      
set(gca,'xscale','log','yscale','log')
title(title_nfib)
xlabel('\Phi')
ylabel('\sigma_{0} L^4/E_YI')
set(gca,'fontsize',14,'fontName','Times New Roman')

sig0_phi_log = log(sig0_phi(1:3));
sig0_log = log(sig0(1:3));
corr = polyfit(sig0_phi_log,sig0_log,1);
sigxz_regressX = linspace(sig0_phi(1), sig0_phi(3),150);
alpha = exp(corr(2));
beta = corr(1);
sigxz_regress = alpha*sigxz_regressX.^beta; 
plot(sigxz_regressX,sigxz_regress,'linewidth',1.3,'color',rgb('MediumBlue'));

else
    
figure(2)
subplot(2,2,loc)
hold on;
box on;
scatter(sig0_phi,sig0,40,'MarkerEdgeColor',rgb('MediumBlue'),...
              'MarkerFaceColor',rgb('MediumBlue'));      
set(gca,'xscale','log','yscale','log')
title(title_nfib)
xlabel('\Phi')
ylabel('\sigma_{0} L^4/E_YI')
set(gca,'fontsize',14,'fontName','Times New Roman')

sig0_phi_log = log(sig0_phi);
sig0_log = log(sig0);
corr = polyfit(sig0_phi_log,sig0_log,1);
sigxz_regressX = linspace(sig0_phi(1), sig0_phi(4),200);
alpha = exp(corr(2));
beta = corr(1);
sigxz_regress = alpha*sigxz_regressX.^beta; 
plot(sigxz_regressX,sigxz_regress,'linewidth',1.3,'color',rgb('MediumBlue'));
    
    
end

end