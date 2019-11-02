lr=0:0.01:4;
L=(1+2*lr).^2./(1+3*lr+3*lr.^2);

figure(1)
plot(lr,L)
xlabel('\itl_r');
ylabel('\itL');
set(gca,'XTickMode','manual','XTick',[0:0.5:4]);
set(gca,'YTick',[1:0.05:1.4]);
grid on

