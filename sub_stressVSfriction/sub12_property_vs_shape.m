%%%
%%% Plot property across shapes
%%%

clc;
clear;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define all relevant cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cases = {'theta0','theta1','theta3','theta6','helical'};
titleName = {'$\theta_{eq} = 0$','$\theta_{eq} = 0.1$',...
    '$\theta_{eq} = 0.3$','$\theta_{eq} = 0.6$',...
    '$(\theta_{eq},\phi_{eq})=(0.8,0.7)$'}; 
nCases = length(cases); 


nfibArr = [160 240 320 640 1280 3200 6400];
nNfib = length(nfibArr); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load all relevant cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nCases
    name = (['sub1_',cases{i},'.mat']); 
    C1(i) = load(name); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up plotting for each Nfib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine legend
legendNfib = cell(nNfib,nCases); 
caseplotted = zeros(nNfib,1); 
for j=1:nNfib
    ind = 1; 
    for i=1:nCases
        for jj=1:C1(i).nNfib
            if C1(i).nfibArr(jj) == nfibArr(j)
                legendNfib{j,ind} = titleName{i};
                ind = ind + 1; 
                break; 
            end
        end
    end
    caseplotted(j) = ind-1; 
end

for jj=1:nNfib
    figure('units','normalized','outerposition',[0.1 0.1 0.6 0.6])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I vs. shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nCases
    for j=1:C1(i).nNfib
        % determine which figure to put it
        for jj=1:nNfib
            if C1(i).nfibArr(j) == nfibArr(jj)
                figure(jj)
                hold on
                muArr = C1(i).muArr;
                A = C1(i).Eelas(:,j);
                plot(muArr,A,'-.o','MarkerSize',10,'Linewidth',2.5);
                break
            end
        end
    end
end

for j=1:nNfib
    figure(j)
    hold on
    title(['$N_{fib} =$ ',num2str(nfibArr(j))])
    xlabel('$\mu$')
    %     ylabel('$\bar{I}$')
    ylabel('$<E_{elastic}>/ 8\pi\eta_0\dot{\gamma}l^3$')
    legend(legendNfib{j,1:caseplotted(j)},'location','bestoutside')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sigP vs. shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:nCases
%     for j=1:C1(i).nNfib
%         % determine which figure to put it
%         for jj=1:nNfib
%             if C1(i).nfibArr(j) == nfibArr(jj)
%                 figure(jj)
%                 hold on
%                 muArr = C1(i).muArr;
%                 particle_stress = C1(i).particle_stress(:,j);
%                 scatter(muArr,particle_stress,'filled');
%                 break
%             end
%         end
%     end
% end
% 
% for j=1:nNfib
%     figure(j)
%     hold on
%     title(['$N_{fib} =$ ',num2str(nfibArr(j))])
%     xlabel('$\mu$')
%     ylabel('$\bar{\sigma_{xz}}$')
%     legend(legendNfib{j,1:caseplotted(j)},'location','bestoutside') 
% end
