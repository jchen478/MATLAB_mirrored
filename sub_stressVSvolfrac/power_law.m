%% load relevant simulation parameters
simulation_cases;

data = load('../sub_redispersion/Basis.mat');
etaData = data.etaData;
NCData = data.NCData;


%% Find power law parameter
%  For every friction and attraction combination
power = zeros(nMu,nAtt,2);
weightFracLog = log(weightfracArr.value)';
for i=1:nMu
    for k=1:nAtt
        display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k))])
        etalog = squeeze(log(etaData(i,k,:)));
        if isinf(min(etalog))
            display('Contains Inf')
            etalog
            continue;
        end
        % find slope and R2 interval
        [power(i,k,2), power(i,k,1)] = regression(weightFracLog(1:end)',etalog(1:end)');
        
    end
end

figure('Units','Inches','Position',[1 1 3.5 3.0]);
box on
hold on;
for k=1:nAtt-2
    plot(muC.value,power(:,k),'-o')
end
legend(attC.legend,'location','best')
xlabel(muC.name)
ylabel('$n$')
