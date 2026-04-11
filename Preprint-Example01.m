y = linspace(0.1,100,100);
z = 3.*log(y) - 3;

figure(1)
plot(y,z)
title('Plot of Restricted Log-Likelihood Function','fontsize',18)
xlabel('\beta','interpreter','tex','fontsize',14)
ylabel('L(1, \beta; (1, 1, 1)^{T})','interpreter','tex','fontsize',14)
