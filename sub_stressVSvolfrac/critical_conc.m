%% finid critical concentrations as function of rp and n (number of contacts)
close all
rp = caseArr('$r_p$',65:10:85,1);
n = caseArr('$n$',1.8:0.1:4,1);

figure('Units','Inches','Position',[1 1 3.5 2.5]);
hold on
for i=1:rp.ndata
    cv = 2*pi/rp.value(i)^2*n.value.^3./(n.value-1);
    plot(cv*100,n.value);
end
legend(rp.legend,'location','best')
xlabel('$c_v (\%) $')
ylabel('$n$')
set(gca,'XMinorTick','on')