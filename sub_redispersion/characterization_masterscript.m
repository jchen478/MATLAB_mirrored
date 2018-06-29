%% analysis for redispersion - include intensity and cluster analysis
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 0 - definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%
simulation_cases;
basisStrain = 1000;
IData = zeros(nMu,nAtt,nA);
IDataB = zeros(nMu,nAtt);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part I - obtain basis
%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 0; 
for i=1:nMu
    for k=1:nAtt
        
        display(['Processing basis mu',num2str(muArr(i)),'_att',num2str(attArr(k))]);
        filePrefix = [dataPath,shape,'_basis_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_'];

        %% read input and process data
        IDataB(i,k) = process_intensity(filePrefix, basisStrain, 2, r);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part II - obtain redispersed value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nMu
    for k=1:nAtt
        for j=1:nA
            display(['Processing mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j))])
            filePrefix = [dataPath,shape,'_mu',num2str(muArr(i)),'_att',num2str(attArr(k)),'_a',num2str(aArr(j)),'_'];
            
            %% read input and process data
            [box_strain, sidex] = read_box([filePrefix,'Lbox.txt']);
            % define intervals
            r = round(intervals(box_strain,sidex),0)';
            % intensity energy statistics
            IData(i,k,j) = process_intensity(filePrefix, basisStrain, 5, r);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part III - calculate property difference due to redispersion cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DIData = zeros(nMu,nAtt,nA);
IData (IData == 0) = NaN;

for i=1:nA
    DIData(:,:,i) = IData(:,:,i) - IDataB;
end

IDataObj = data3D('$I$',IData); 
DIDataObj = data3D('$\Delta I$',DIData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intensity
plot3dim1(IDataObj,attC,muC,aC); 
plot3dim2(DIDataObj,IDataObj,attC,muC,aC)
