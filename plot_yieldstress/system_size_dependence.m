clc;
close all;

	%% beta			
beta = [320 4.824   0.407	
    640 4.119 0.116
    1280	4.958	0.098	;
	2560	4.493	0.221	;
	6400	4.170	0.198	;
	12800	4.753	0.301	];
				
	%% alpha			
alpha = [320 3.20E10 6.53E00	
    640 2.11E9 1.453
    1280	1.29E+11	1.32E+00	;
	2560	2.06E+10	2.36E+00	;
	6400	8.37E+09	2.22E+00	;
	12800	1.14E+11	3.50E+00	];

figure(1)
hold on;
errorbar(beta(:,1),beta(:,2),beta(:,3)/2,'o','color',rgb('MediumBlue'))
xlabel('N_{fib}')
ylabel('\beta')
set(gca,'fontsize',22,'fontName','Times New Roman')

figure(2)
hold on;
scatter(alpha(:,1),alpha(:,2),40,'MarkerEdgeColor',rgb('MediumBlue'),...
              'MarkerFaceColor',rgb('MediumBlue'))
xlabel('N_{fib}')
ylabel('\alpha (Pa)')
set(gca,'fontsize',22,'fontName','Times New Roman')
