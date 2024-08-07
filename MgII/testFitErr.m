clear
x = linspace(0,10,8);
xe = randn([1,8])/10;
ye = randn([1,8])/10;
y = 6*(x+xe) +10 +ye;
yee = 6*xe + ye;

[a_coef_dist, b_coef_dist, R, a_coef, b_coef] = errFitter(x,y,xe,yee, 1000);

fig = figure();
errorbar(x,y,yee, 'LineStyle','None')
hold on 

for i=1:1000
    p=plot(x, a_coef_dist(i)*x + b_coef_dist(i));
    p.Color = [1,0,1, 0.02];
    p.LineWidth= 0.1;
    hold on
end
title(sprintf('%.3f_{%.3f}^{%.3f}x + (%.3f)_{%.3f}^{%.3f}, <R>=%.3f', a_coef, b_coef, median(R)))
xlabel('x')
ylabel('y = 6x+10','interpreter','latex')
exportgraphics(fig, 'testFiterr.png', 'Resolution',800)