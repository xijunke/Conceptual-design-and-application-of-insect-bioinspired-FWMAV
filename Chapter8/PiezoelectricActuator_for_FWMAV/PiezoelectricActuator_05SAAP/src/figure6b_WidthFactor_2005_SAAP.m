% w_r=0:0.01:2;
w_r=0.01:0.02:1.99;
%  A=8*(1-w_r).^3;
%  B=3*log((2-w_r)./w_r).*(w_r-2).^2;
%  C=-12*w_r.^2+30*w_r-18;
%  W=A./(B+C);
 W=8*(1-w_r).^3./(3*(w_r-2).^2.*log((2-w_r)./w_r)-18+30*w_r-12*w_r.^2);
 
 figure(2)
 plot(w_r,W)
 xlabel('\itw_r');
ylabel('\itW');
set(gca,'XTickMode','manual','XTick',[0:0.25:1.99]);
set(gca,'YTick',[0:0.1:1.4]);
grid on





