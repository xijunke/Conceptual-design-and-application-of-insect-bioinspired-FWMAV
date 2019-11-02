clear all;clc;
x=0.01:0.01:1.99;
y=0:0.1:4;
[wr lr]=meshgrid(x,y);
gc=8*(1-wr).^3;
gd=6*(wr-1).*(3+4*lr-2*wr-4*lr.*wr);
ge=3*(-2-2*lr+wr+2*lr.*wr).^2.*log((2-wr)./wr);
U=(1+2*lr).^2.*gc./(gd+ge);

figure(1)
% subplot(121)
mesh(wr,lr,U);
xlabel('\itw_r');
ylabel('\itl_r');
zlabel('\itG_U');

figure(2)
% subplot(122)
[C,h]=contour(wr,lr,U);
clabel(C,h)
xlabel('\itw_r');
ylabel('\itl_r');