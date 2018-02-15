load('nfibArr.mat')

figure('units','normalized','outerposition',[0.1 0.1 0.5 0.6])
histogram(nfibArr,10,'EdgeColor',rgb('MidnightBlue'),'FaceColor',rgb('DarkBlue'))
xlabel('$N_{fib}$')
title('Number of fibers used in Switzer''s simulations')

print('nfibdistribution','-dpng')