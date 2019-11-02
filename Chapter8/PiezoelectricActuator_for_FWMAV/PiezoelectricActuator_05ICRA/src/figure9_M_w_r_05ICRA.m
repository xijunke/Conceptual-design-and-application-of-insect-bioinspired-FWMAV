clear all;
clc;
M_wr=[];
for wr=0:0.001:2;
%wr=0.5;
g1=wr-1;
g2=wr-2;
g4=2*g1*(-39624+wr*(135808+wr*(-182782+wr*(120878+wr*(-39257+4992*wr)))));
g5=450*g2^6*log(2-wr)^2;
g6=15*g2^4*log(2-wr)*(-396+(676-279*wr)*wr-60*g2^2*log(wr));
g7=15*g2^4*log(wr)*(396+(-676+279*wr)*wr+30*g2^2*log(wr));
g8=3600*g1*(-6+10*wr-4*wr^2+g2^2*log((2-wr)/wr));
mw=(g4+g5+g6+g7)./g8;
M_wr=[M_wr mw];
end
wr=0:0.001:2;
figure(1)
plot(wr,M_wr)
xlabel('\itw_r');
ylabel('\itM_{wr}')
% axis([0,2,0.05,0.5])
