t=0:0.1:90;
y=t*pi/180;
Q=atan(0.0694*sin(2*y));
T=Q*180/pi;
plot(t,T)