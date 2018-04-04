%%%
%%% Plot transient percolation results
%%%

clc;
clear;
close all;

% Define common parameters
simulation_cases;

% data path
dataPath = '../data_stressVSfriction/ClusterP/';

% replicate info
nReplicate = 2;
replicate_flag = zeros(nNfib,1);

% figure info
figStart = 1;
markersize = 10;
dataFile = cell(2,1);
for ii=1:2
    dataFile{ii} = ['fig16_principal_replicates/',fileNameArr{1},'_nfib'];
    for j=1:nNfib
        dataFile{ii} = [ dataFile{ii},num2str(nfibArr(j)),'_'];
    end
end
dataFile{1} = [ dataFile{1},'principal_direction3'];

% sampling info
sample_strain = 1000:25:1500;
nSample = length(sample_strain);

% array allocation
Np_avg = zeros(nMu, nNfib);
Np_std = zeros(nMu, nNfib);

%%  Plot principal directions vs. strain

for j=1:nNfib
    figure('Units','Inches','Position',[3 4 3.5 3.5])
    hold on
    title([fileNameArr{1}, ' ', num2str(nfibArr(j))])
    N = zeros(nMu, nReplicate+1);
    for i=1:nMu
        name = [dataPath,fileNameArr{1},'_clusterp_nfib',num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
        File=fopen(name,'r');
        data = fscanf(File,'%f',[36 Inf])';
        fclose(File);
        
        % eliminate non-percolating structures
        sidexo = data(:,5);
        mle3 = data(:,35);
        keep = (mle3 >= sqrt(3)*sidexo);
        
        % delete empty rows
        data = data.*keep;
        data( ~any(data,2),:) = [];
        
        % define variables
        nfib = data(:,1);
        mu = data(:,2);
        strain = data(:,3);
        fcut = data(:,4);
        sidexo = data(:,5);
        sidexn = data(:,6);
        dx = data(:,7);
        nx = data(:,8);
        cID = data(:,9);
        nC = data(:,10);
        rcmx = data(:,11);
        rcmy = data(:,12);
        rcmz = data(:,13);
        Rxx = data(:,14);
        Rxy = data(:,15);
        Rxz = data(:,16);
        Ryy = data(:,17);
        Ryz = data(:,18);
        Rzz = data(:,19);
        x1 = data(:,20);
        y1 = data(:,21);
        z1 = data(:,22);
        x2 = data(:,23);
        y2 = data(:,24);
        z2 = data(:,25);
        x3 = data(:,26);
        y3 = data(:,27);
        z3 = data(:,28);
        lam1 = data(:,29);
        lam2 = data(:,30);
        lam3 = data(:,31);
        Rg = data(:,32);
        mle1 = data(:,33);
        mle2 = data(:,34);
        mle3 = data(:,35);
        % calculate invariants of gyration tensor
        I1 = (lam1+lam2+lam3)./lam3;
        I2 = (lam1.*lam2+lam2.*lam3+lam1.*lam3)./(lam3.^2);
        I3 = (lam1.*lam2.*lam3)./lam3.^3;
        
        % shape measurements
        b = lam3.^2 - 0.5*(lam2.^2+lam1.^2); %asphericity
        c = lam2.^2 - lam1.^2; % acylindricity
        k2 = 3/2*(lam1.^4+lam2.^4+lam3.^4)./(lam1.^2+lam2.^2+lam3.^2).^2-0.5;
        % relative shape anisotropy
        
        % obtain number of percolating structures
        allN = zeros(nSample,1);
        for ii=1:nSample
            tstrain = sample_strain(ii);
            allN(ii) = sum(strain == tstrain);
        end
        N(i,1) = mean(allN);
        
        for r=1:nReplicate
            replicate_name = ['../data_stressVSfriction/ClusterP_replicate',...
                num2str(r),'/',fileNameArr{1},'_clusterp_nfib',...
                num2str(nfibArr(j)),'_',num2str(muArr(i)),'.txt'];
            
            if exist(replicate_name, 'file') == 2
                
                % if replicates exist, open file ...
                replicate_flag(j) = r;
                File=fopen(replicate_name,'r');
                data = fscanf(File,'%f',[35 Inf])';
                fclose(File);
                
                % eliminate non-percolating structures ...
                sidexo = data(:,5);
                mle3 = data(:,35);
                keep = (mle3 >= sqrt(3)*sidexo);
                
                % delete empty rows ...
                data = data.*keep;
                data( ~any(data,2),:) = [];
                
                % obtain strain 
                strain = data(:,3);
                
                % obtain percolation results
                allN = zeros(nSample,1);
                for ii=1:nSample
                    tstrain = sample_strain(ii);
                    allN(ii) = sum(strain == tstrain);
                end
                N(i,1+r) = mean(allN);
            end
        end
        
        % different plots
        % 1. Third principal axis
        %{
        scatter3(abs(x3),abs(y3),abs(z3),'filled')
        box on
        xlabel('$\nu_{3,x}$ (flow)')
        ylabel('$\nu_{3,y}$ (vorticity)')
        zlabel('$\nu_{3,z}$ (gradient)')
        %}
        % 2. mle3 * Third principal axis
        %{
        scatter3(mle3.*abs(x3),mle3.*abs(y3),mle3.*abs(z3),'filled')
        box on
        xlabel('$\nu_{3,x}MLE_3$ (flow)')
        ylabel('$\nu_{3,y}MLE_3$ (vorticity)')
        zlabel('$\nu_{3,z}MLE_3$ (gradient)')
        %}
        % 3. mle3
        %{
        scatter(mu,mle3,'filled')
        xlabel('$\mu$')
        ylabel('$MLE_3$')
        %}
        % 4. mle2
        %{
        scatter(mu,mle2,'filled')
        xlabel('$\mu$')
        ylabel('$MLE_2$')
        %}
        % 5. Rg
        %{
        scatter(mu,Rg,'filled')
        xlabel('$\mu$')
        ylabel('$R_g$')
        %}
        % 6. I1
        %{
        scatter(mu,I1,'filled')
        xlabel('$\mu$')
        ylabel('$I_1$')
        %}
        % 7. I2
        %{
        scatter(mu,I2,'filled')
        xlabel('$\mu$')
        ylabel('$I_2$')
        %}
        % 8. I3
        %{
        scatter(mu,I3,'filled')
        xlabel('$\mu$')
        ylabel('$I_3$')
        %}
        % 9. strain vs. data
        %{
        scatter(strain,c./(3*sidexo).^2,markersize,...
            'MarkerFaceColor',colorArr{i},...
            'MarkerEdgeColor',colorArr{i})
        xlabel('$\gamma$')
        ylabel('$c/(3L_{box})^2$')
        xlim([1000 1500])
        %}
        scatter(strain,nC,markersize,...
            'MarkerFaceColor',colorArr{i},...
            'MarkerEdgeColor',colorArr{i})
        xlabel('$\gamma$')
        %         ylabel('$mle1$ in floc')
        
    end
    % average over replicates
    average_range = replicate_flag(j) + 1;
    Np_avg(:,j) = mean(N(:,1:average_range),2);
    Np_std(:,j) = std(N(:,1:average_range),0,2);
end

%% plot number of percolation structures
figure('Units','Inches','Position',[0.5 0.5 6.5 3.5])
ylabel('$<N_p>$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,Np_avg(:,j),...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    else
        errorbar(muArr,Np_avg(:,j),...
            Np_std(:,j)/2,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end
xlabel('$\mu$')
xlim([0 inf])

legend(thetaNfibLegendArr,'location','bestoutside')

figure('Units','Inches','Position',[0.5 4 6.5 3.5])
ylabel('$<N_p> L/L_{box}$')
hold on
for j=1:nNfib
    if replicate_flag(j) == 0
        plot(muArr,Np_avg(:,j)/lboxArr(j)*2*rpFiber,...
            '-.o','MarkerSize',10, 'linewidth',2.5);
    else
        errorbar(muArr,Np_avg(:,j)/lboxArr(j)*2*rpFiber,...
            Np_std(:,j)/2/lboxArr(j)*2*rpFiber,...
            '-.o','MarkerSize',8, 'linewidth',2.5);
    end
end
xlabel('$\mu$')
xlim([0 inf])
legend(thetaNfibLegendArr,'location','bestoutside')
