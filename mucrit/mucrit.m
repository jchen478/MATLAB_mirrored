clc;
% clear all;
close all;
%% parameters
dx = [29.75 59.50 79.33 119.00 238.00];
kb = [6 8 10 12];
mu = [2.0 2.5 3.0 3.5 4.0 4.5 5.0];

n_dx = length(dx);
n_kb = length(kb);
n_mu = length(mu);
n_cases = n_dx*n_kb*n_mu;

window = 200;
time = 0:0.2:1500;
n_data = length(time); 
Hts = zeros(n_data,n_cases); 
H = zeros(n_data,n_cases);
F = cell(n_cases,1);
S = cell(n_cases,1);

%% read data
for i=1:n_dx
    for j=1:n_kb
        for k=1:n_mu
            
            ind = (i-1)*n_kb*n_mu+(j-1)*n_mu+k;
            F{ind} = ['Shannon_Entropy_',num2str(dx(i),'%.2f'),'_',num2str(kb(j)),'_',num2str(mu(k),'%.1f'),'.txt'];
            S{ind} = ['dx = ',num2str(dx(i),'%.2f'),', kb = ',num2str(kb(j)),', \mu = ',num2str(mu(k),'%.1f')];
            File = fopen(F{ind},'r');
            data = fscanf(File,'%f %f',[2 inf])';
            fclose(File);
            H(:,ind) = data(2:end,2);
            
        end
    end
end

%% time-smooth data
for i=1:n_cases
    
    for j=1:n_data-window
       
        sum = 0;
        for k=1:window
           
            sum = sum + H(j+k,i);
            
        end
        sum = sum/window; 
        Hts(j,i) = sum;
    end
    for j=n_data-window+1:n_data
        Hts(j,i) = mean(H(j:end,i)); 
    end
end

%% plot

% one plot for every dx
% for i=1:n_dx
%     figure(i)
%     hold on;
%     for j=1:n_kb     
%         for k=1:n_mu
%             ind = (i-1)*n_kb*n_mu+(j-1)*n_mu+k;
%             plot(time,smoothts(H(:,ind),'b',2000));
%         end  
%     end
%     llim = (i-1)*n_kb*n_mu+1;
%     hlim = i*n_kb*n_mu ;
%     legend(S{llim:hlim});
%     title(['dx = ',num2str(dx(i))])
%     set(gcf,'color','white')
%     set(gca,'fontsize',24,'fontName','Times New Roman')
% end

% one plot per kb for selected dx
i = 4;
for j=1:n_kb
    figure(j)
    hold on;
    for k=1:n_mu
        ind = (i-1)*n_kb*n_mu+(j-1)*n_mu+k;
        plot(time,Hts(:,ind));
    end
    llim = (i-1)*n_kb*n_mu+(j-1)*n_mu+1;
    hlim = (i-1)*n_kb*n_mu+(j)*n_mu;
    legend(S{llim:hlim});
    title(['dx = ', num2str(dx(i)),', kb = ',num2str(kb(j))])
    set(gcf,'color','white')
    set(gca,'fontsize',24,'fontName','Times New Roman')
    ylim([3.95 4.15])
end




