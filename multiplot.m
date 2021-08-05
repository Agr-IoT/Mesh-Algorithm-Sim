graphics_toolkit("gnuplot") 

x = linspace(0,10,50);
y1 = sin(x);
plot(x,y1)
title('Combine Plots')

hold on

y2 = sin(x/2);
plot(x,y2)

y3 = 2*sin(x);
scatter(x,y3) 

hold off