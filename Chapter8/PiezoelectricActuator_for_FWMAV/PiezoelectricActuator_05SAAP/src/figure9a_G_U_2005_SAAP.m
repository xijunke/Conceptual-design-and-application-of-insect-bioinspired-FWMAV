 x=0.04:0.05:1.99;
 y=0.09:0.1:3.99;
 [w_r,l_r]=meshgrid(x,y);
 g_delta=(1+2*l_r);
 g_c=8*(1-w_r).^3;
 g_d=6*(w_r-1).*(3+4*l_r-2*w_r-4*l_r.*w_r);
 g_e=3*(-2-2*l_r+w_r+2*l_r.*w_r).^2.*log((2-w_r)./w_r);
 G_le=(g_d+g_e)./g_c;
 G_F=g_delta./G_le;
 G_U=g_delta.*G_F;
 
figure(1)
surf(w_r,l_r,G_U)
xlabel('\itw_r');
ylabel('\itl_r');
zlabel('\itG_U');