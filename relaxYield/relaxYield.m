function  [SteadyValPreRelax,SteadyVal,SeffArr] = relaxYield(File_Parameters,File_Stress,File_NC,File_Eelastic,figNum,strainInd,SteadyValPreRelax,SteadyVal,startAvg,SeffArr)

%% Parameters
File = fopen(File_Parameters,'r');
Parameters = textscan(File,'%f %*s',Inf,'Delimiter','\n');
fclose(File);
Parameters = cell2mat(Parameters);

%% Stress tensor
File = fopen(File_Stress,'r');
Stress = fscanf(File,'%f',[7 Inf])';
fclose(File);

%% Number of Contacts
File = fopen(File_NC,'r');
NC = fscanf(File,'%f',[5 Inf])';
fclose(File);

%% Elastic storage
File = fopen(File_Eelastic,'r');
elastic = fscanf(File,'%f',[4 Inf])';
fclose(File);
elastic(1,:) = []; 

%% Define parameters
nfib = Parameters(1); 
nseg = Parameters(2); 
rps  = Parameters(3); 
kb  = Parameters(4); 
dt = Parameters(10);
strain = Parameters(11); 
side = Parameters(12); 
config_write = Parameters(14); 
contact_write = Parameters(15); 
stress_write = Parameters(49); 
strain_stop = Parameters(68); 
elastic_write = Parameters(69);

%% Calculate parameters
L = 2*rps*nseg; 
nL3 = nfib*L^3/side^3;
Seff = kb*pi/nseg^4;
over_Seff = 1/Seff; 

SeffArr(figNum) = Seff; 

%% Time 
tStress = dt*stress_write:dt*stress_write:strain;
tCon = dt*contact_write:dt*contact_write:strain;
telastic =  dt*elastic_write:dt*elastic_write:strain;

tStressStop = tStress-strain_stop;
tConStop = tCon-strain_stop;
telasticStop = telastic-strain_stop; 

%% Stress processing
sigxzStar = Stress(:,4)*pi*nL3/(6*nseg^3*log(2*rps))*over_Seff + over_Seff;

%% Figure formatting parameters
markersize = 10;

%% Plotting
figure(figNum)
subplot(3,1,1)
hold on;
scatter(tStressStop,sigxzStar,markersize)
ylabel('\it{\sigma_{xz} L^4/ E_Y I}')

subplot(3,1,2)
hold on;
scatter(tConStop,NC(:,4),markersize)
ylabel('N_C / N_{fib}')

subplot(3,1,3)
hold on;
    scatter(telasticStop,elastic(:,4),markersize)
ylabel('\it{}E_{el}*')

%% calculate steady state values before relaxation 
nendStress = strain_stop/stress_write/dt;
nendCon = strain_stop/contact_write/dt;
nendElastic = strain_stop/elastic_write/dt;
SteadyValPreRelax(strainInd,figNum,1) = mean(sigxzStar(1:nendStress));
SteadyValPreRelax(strainInd,figNum,2) = mean(NC(1:nendCon,4));
SteadyValPreRelax(strainInd,figNum,3) = mean(elastic(1:nendElastic,4));

%% calculate steady state values
start = strain_stop + startAvg;
nstartStress = start/stress_write/dt;
nstartCon = start/contact_write/dt;
nstartElastic = start/elastic_write/dt;
SteadyVal(strainInd,figNum,1) = mean(sigxzStar(nstartStress:end));
SteadyVal(strainInd,figNum,2) = mean(NC(nstartCon:end,4));
SteadyVal(strainInd,figNum,3) = mean(elastic(nstartElastic:end,4));

end
