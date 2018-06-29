n = logspace(1,5);

d = n.^0.2; 
e = (log10(n)).^5; 

plot(n,d)
hold on
plot(n,e)

legend('d','e')