function sizeCharFun(filename,markersize,linewidth,color,floc_cutoff)    


%% Open and read file
File = fopen(filename,'r');
dataLarge = fscanf(File,'%f',[18 Inf])';
fclose(File);

%% Allocate space for flocs exceeding cutoff
nfibFloc = dataLarge(:,14);
keep = nfibFloc >= floc_cutoff;
data = zeros(sum(keep),18);

%% Eliminate rows with number of fibers below cutoff
ind = 1; 
for i=1:size(dataLarge,1)
    if (keep(i) == 1)
        data(ind,:) = dataLarge(i,:); 
        ind = ind + 1;
    end
end


%% Assign variables
nfib = data(:,1);
nseg = data(:,2);
mustat = data(:,3);
mukin = data(:,4);
rps = data(:,5);
sidex = data(:,6);
sidey = data(:,7);
sidez = data(:,8);
dx = data(:,9);
dy = data(:,10);
dz = data(:,11);
flocId = data(:,12);
% floc_cutoff = data(:,13);
nfibFloc = data(:,14);
maxGap = data(:,15);
sizex = data(:,16);
sizey = data(:,17);
sizez = data(:,18);

%% mean cluster size calculation
S = sum(nfibFloc.^2)/sum(nfibFloc)
Sx = sum(sizex.^2)/sum(sizex)

%% size
figure(1)
subplot(1,3,1)
hold on
scatter(mustat,sizex,markersize,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'Linewidth',linewidth);
ylabel('$S_x$')
% title('Size in x')

subplot(1,3,2)
hold on
scatter(mustat,sizey,markersize,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'Linewidth',linewidth);
ylabel('$S_y$')
% title('Size in y')
% title(['Size of floc ( Max Gap = ',num2str(maxGap(1)), ',   Floc cutoff = ',num2str(floc_cutoff),')'])
title(['Size of floc ( Box size: ',num2str(sidex(1)),' x ',num2str(sidey(1)), ' x ', num2str(sidez(1)),', Floc cutoff = ',num2str(floc_cutoff),')'])

subplot(1,3,3)
hold on
scatter(mustat,sizez,markersize,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'Linewidth',linewidth);
ylabel('$S_z$')
% title('Size in z')

%% fraction
figure(2)
subplot(1,3,1)
hold on
scatter(mustat,sizex./sidex,markersize,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'Linewidth',linewidth);
ylabel('$S_x / L_x$')
% title('Fraction in x')

subplot(1,3,2)
hold on
scatter(mustat,sizey./sidey,markersize,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'Linewidth',linewidth);
ylabel('$S_y / L_y$')
% title('Fraction in y')
title(['Size of floc / box length ( Max Gap = ',num2str(maxGap(1)), ',   Floc cutoff = ',num2str(floc_cutoff),')'])

subplot(1,3,3)
hold on
scatter(mustat,sizez./sidez,markersize,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'Linewidth',linewidth);
ylabel('$S_z / L_z$')
% title('Fraction in z')

figure(3)
subplot(1,3,1)
hold on
scatter(nfibFloc,sizex,markersize,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'Linewidth',linewidth);
xlabel('Number of fibers in floc')
ylabel('$S_x$')

subplot(1,3,2)
hold on
scatter(nfibFloc,sizey,markersize,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'Linewidth',linewidth);
xlabel('Number of fibers in floc')
ylabel('$S_y$')
% title(['Size of floc ( Max Gap = ',num2str(maxGap(1)), ',   Floc cutoff = ',num2str(floc_cutoff),')'])

subplot(1,3,3)
hold on
scatter(nfibFloc,sizez,markersize,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'Linewidth',linewidth);
xlabel('Number of fibers in floc')
ylabel('$S_z$')

end